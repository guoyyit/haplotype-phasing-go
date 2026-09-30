# 输入/输出文件格式

所有文件均为 UTF-8、制表符分隔（TSV）。空行与以 `#` 开头的行会被忽略；
单独成行的 `#variants` / `#genotypes` / `#reads` / `#trios` 是段标记，
可用于自描述。列标题行使用下划线形式的标识（如 `variant_id`、`sample_id`、
`read_id`、`child_id`），会被自动识别并跳过。

位点按文件中的顺序赋予 0 基下标，分相、块与输出均按该顺序。

## 1. 位点文件（`-variants`）

每行：`variant_id chrom pos ref alt`

| 列 | 含义 |
|---|---|
| variant_id | 全局唯一位点 ID |
| chrom | 染色体 |
| pos | 非负整数坐标 |
| ref | 参考等位基因（规范 A0） |
| alt | 备选等位基因（规范 A1） |

仅支持二等位基因位点。示例见 `examples/variants.tsv`。

## 2. 基因型文件（`-genotypes`）

每行：`sample variant_id gt`

- 样本按首次出现顺序记录，输出沿用该顺序。
- `gt` 为未分相基因型，分隔符 `/` 或 `|` 等价（输入的相位不被采信）。
- 缺失与部分缺失：
  - `./.`：两条等位基因均缺失（倍性 0）。
  - `A/.`：仅观测到一条等位基因（倍性 1）。
  - `A/A`：纯合；`A/G`：杂合。
- 未出现在文件中的“样本×位点”组合按完全缺失处理，不会导致失败。

规范等位基因顺序：杂合位点的两个不同等位基因以 `ref` 为 `A0`、`alt` 为
`A1`。输出中的 `h1`/`h2` 用 `A0`/`A1` 表示该单倍型携带哪个规范等位基因。

## 3. 读段文件（`-reads`，可选）

每行：`read_id sample variant_id allele`

- 相同 `(read_id, sample)` 的多行按位点顺序聚合成一条读段。
- `allele` 必须等于该位点的 `ref` 或 `alt`；`.` 为空观测。
- 无效观测（未知位点、缺失/非 ref-alt 等位基因）被忽略，打印 `notice`
  并计入冲突/告警，不会静默丢弃或中止。
- 同一条读段在同一位点出现多个矛盾等位基因时，按等位基因字典序取最小者，
  该“读段内部冲突”会作为 notice 与冲突计数上报。

## 4. 家系文件（`-trios`，可选）

每行：`child father mother`

每个 child 至多一行。父、母、子必须是基因型文件中的样本 ID；缺失样本不
产生遗传约束（视为信息不足），不报错。

## 5. 分相输出文件（`-out`，格式版本 `1.0`）

段按固定顺序排列，便于机读与重读：

```
#haplophased  1.0
#variants        位点定义（同输入）
#trios           家系定义（同输入）
#reads           全部有效读段观测（保证可复现读段证据）
#config          单行紧凑 JSON，即本次运行使用的完整配置
#phase           每样本每位点一行的分相结果
#blocks          每样本每个分相块一行
#stats           运行级统计键值
```

`#phase` 列：

| 列 | 含义 |
|---|---|
| sample | 样本 ID |
| variant_index | 0 基位点下标 |
| variant_id | 位点 ID |
| status | `PHASED` / `UNPHASED` / `HOMOZYGOUS` / `PARTIAL` / `MISSING` |
| low_cov | `1`=覆盖度低于 `min_site_depth` |
| a0, a1 | 规范等位基因（ref/alt） |
| h1, h2 | 两条单倍型的等位基因：`A0`、`A1` 或 `.`（未分相/缺失） |
| violation | `1`=违反家系孟德尔遗传约束 |

status 语义：

- `PHASED`：杂合位点且已分相，`h1`/`h2` 为两个不同等位基因。
- `UNPHASED`：杂合位点但未分相，或违反家系约束（同时 `violation=1`），
  不计入正确分相。
- `HOMOZYGOUS`：纯合位点，两条单倍型相同。
- `PARTIAL`：仅观测到一条等位基因。
- `MISSING`：完全缺失基因型。

`#blocks` 列：`sample block_id start_idx end_idx site_count hets_count
support conflict confidence error_rate status`。

- `block_id` 形如 `BLK-<sample>-0001`，按样本与位置确定性编号。
- `start_idx`/`end_idx` 为位点下标闭区间。
- `support`/`conflict` 为块内已裁决边的多数/少数票合计。
- `confidence` 与 `error_rate` 均四舍五入到 4 位小数且和为 1。
- status：`OK`、`LOW_CONF`（置信度低于 `min_block_confidence`，仍输出）、
  `SINGLETON`（单个杂合位点，无相邻证据，置信度记 1）。

`#stats` 的键与 `model.RunStats` 一一对应，含义见
[algorithm.md](algorithm.md) 的“统计字段”。

## 6. 报告文件（`-report`）

人读纯文本：总体计数与每样本的块、置信度、违规摘要。该文件仅供阅读，
不参与重读；机读与往返校验以 `-out` 文件为准。
