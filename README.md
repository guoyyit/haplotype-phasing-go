# haplotyper

本地运行的遗传变异单倍型分相服务。输入同一批二倍体常染色体二态变异位点上多个样本的
未分相基因型，可选提供跨位点测序读段与父/母/子三人家系，输出每个样本的两条单倍型、
分相块、块置信度与全部证据冲突计数。所有输入输出均为本地文件，无任何网络依赖。

## 构建

```sh
go build ./...
```

需要 Go 1.25+。

## 运行

```sh
go run ./cmd/haplotyper \
  -variants variants.tsv \
  -genotypes genotypes.tsv \
  -reads reads.tsv \          # 可选
  -pedigree trio.tsv \        # 可选
  -config config.json \       # 可选
  -output result.json
```

回读既有结果并重新分相（输出须与原结果逐条一致）：

```sh
go run ./cmd/haplotyper -variants V -genotypes G -reads R -pedigree P \
  -rephase result.json -output result.again.json
diff result.json result.again.json   # 无差异
```

参数：

- `-variants`（必填）：变异位点 TSV。
- `-genotypes`（必填）：样本基因型 TSV。
- `-reads`（可选）：测序读段观测 TSV。
- `-pedigree`（可选）：恰好一个父/母/子三人 TSV。
- `-config`（可选）：JSON 配置，缺省字段使用内置默认值。
- `-output`（必填）：结果 JSON 输出路径。
- `-rephase`（可选）：此前由本服务写出的结果 JSON。

## 输入文件格式

所有 TSV 使用制表符分隔；空行与以 `#` 开头的行被忽略。

### 变异位点 `variants.tsv`

表头固定为 `variant_id chrom pos ref alt`：

```tsv
variant_id	chrom	pos	ref	alt
v1	chr1	100	A	T
v2	chr1	200	C	G
```

位点按 `(chrom, pos, variant_id)` 升序确定为内部规范顺序；输出按该顺序给出。

### 基因型 `genotypes.tsv`

表头第一列为 `sample_id`，其余列为变异 `variant_id`（顺序任意，但必须都存在于变异文件）：

```tsv
sample_id	v1	v2
child	0/1	./.
father	0|0	0/1
```

基因型 token：`0/0`、`0/1`、`1/0`、`1/1`，其中 `0` 表示 ref、`1` 表示 alt。
`/` 与 `|` 等价（输入恒视为未分相）。`.` 表示等位缺失：

- `./.` 或单独的 `.`：完全缺失，输出 `missing`。
- `0/.`、`./1`：部分缺失，保留已知等位，输出 `unphased`。

### 测序读段 `reads.tsv`

表头：`read_id sample_id variant_id allele [weight]`。

```tsv
read_id	sample_id	variant_id	allele
r1	child	v1	0
r1	child	v2	1
```

- 同一 `read_id` 的多行属于同一条物理读段，必须来自同一个样本。
- 一条读段在同一个变异上只能有一个等位；重复行须等位一致，权重相加。
- `allele` 仅允许 `0`/`1`；`weight` 缺省为 `1.0`，必须为正数。
- 一条读段覆盖的两个相邻杂合位点构成一条“同源归属”证据（见算法文档）。

### 家系 `trio.tsv`

表头固定为 `child father mother`，恰好一行数据：

```tsv
child	father	mother
child	father	mother
```

三个样本都必须存在于基因型文件中。

### 配置 `config.json`

字段及默认值见“配置项”一节。

## 输出格式 `result.json`

稳定的 JSON（schema 版本 `haplotyper-result/v1`），字段顺序固定，浮点置信度四舍五入到
12 位小数，保证相同输入与配置下逐字节一致。顶层结构：

- `schema_version`：输出版本标识。
- `config`：实际生效的完整配置（含默认值与随机种子）。
- `variants`：规范顺序的位点列表。
- `samples[]`：每个样本的 `blocks` 与 `sites`。
- `conflicts`：冲突计数与逐条明细。
- `stats`：全局限点与分相统计。

`stats` 字段：`samples`、`variant_sites`、`heterozygous_sites`（完全已知杂合位点）、
`phased_hets`（已分相杂合）、`unphased_hets`（未分相杂合）、`partial_sites`（部分缺失位点）、
`missing_sites`（完全缺失位点）、`low_coverage_sites`、`correctly_phased`
（已分相且非孟德尔违规的杂合位点数）。

每个 `sites[]` 条目：

| 字段 | 含义 |
| --- | --- |
| `variant_idx` / `variant_id` | 规范位点序号与 ID |
| `status` | `phased` / `unphased` / `missing` |
| `h1` / `h2` | 两条单倍型等位；未知为 `-1` |
| `block_id` | 分相块编号，`0` 表示不属于任何块 |
| `paternal_h1` | 家系可用时，`h1` 是否为父源等位 |
| `origin_resolved` | 该位点父/母源是否可由家系判定 |
| `mendel_violation` | 是否违反孟德尔遗传约束 |
| `low_coverage` | 读段覆盖数是否低于 `min_reads_per_site` |
| `read_support` | 该位点观测到的读段数 |

每个 `blocks[]` 条目：

| 字段 | 含义 |
| --- | --- |
| `block_id` | 样本内块编号，从 1 开始 |
| `start_variant_idx` / `end_variant_idx` | 块覆盖的位点范围 |
| `heterozygous_sites` | 块内杂合位点数 |
| `supporting_reads` | 块内各边投票读段总数 |
| `confidence` | 块置信度，等于各相邻边后验置信度之积 |
| `expected_error_rate` | `1 - confidence`，即块预期错误率 |

`conflicts` 结构：

- `read_read_conflicts`：投票给被否决奇偶性的读段条数。
- `read_trio_conflicts`：读段奇偶性与家系约束冲突且被家系否决的边数。
- `read_genotype_conflicts`：读段等位与已知纯合基因型矛盾的观测条数（不参与投票）。
- `trio_mendel_violations`：孩子违反孟德尔遗传的位点数。
- `entries[]`：逐条冲突记录（样本、读段、相邻位点、投票奇偶、保留奇偶、类型）。

结果文件可被 `-rephase` 重新读入；详见 `docs/ALGORITHM.md` 的输出回读语义。

## 进一步阅读

- 算法语义、冲突裁决与块置信度：`docs/ALGORITHM.md`
- 全部配置项与边界规模：`docs/CONFIG.md`
- 端到端示例数据：`examples/`
