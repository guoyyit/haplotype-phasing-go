# haplophase-go

本地运行的遗传变异单倍型分相服务。输入为多个样本在同一批二倍体常染色体
二等位基因变异位点上的未分相基因型，可附带覆盖多个位点的测序读段与
父–母–子三人家系；输出每个样本的两条单倍型、分相块划分、每块置信度/
预期错误率，以及读段冲突与家系违规的计数。

## 构建

```sh
go build ./...
go test ./...
```

生成可执行文件：

```sh
go build -o haplophrase ./cmd/haplophrase
```

## 运行

```sh
go run ./cmd/haplophrase \
  -variants examples/variants.tsv \
  -genotypes examples/genotypes.tsv \
  -reads examples/reads.tsv \
  -trios examples/trios.tsv \
  -config examples/config.json \
  -out phased.tsv \
  -report report.txt
```

参数：

- `-variants`（必填）：位点文件。
- `-genotypes`（必填）：未分相基因型文件。
- `-reads`（可选）：测序读段观测文件。
- `-trios`（可选）：父–母–子家系文件。
- `-config`（可选）：JSON 配置；缺省文件或缺省参数使用内置默认值。
- `-out`：分相结果文件（默认 `phased.tsv`）。
- `-report`：人读报告；传空字符串可关闭（默认 `report.txt`）。
- `-reread`：不重新分相，改为读入既往分相结果、用其中保存的输入与配置
  重新计算并做逐字节一致性校验。

读段与家系文件可缺省；缺省时仅依据基因型与确定性默认规则分相。所有输入
输出均为本地文件，服务不访问网络。

## 确定性

相同输入与相同配置多次运行产生逐字节相同的输出。读段处理顺序、边投票
聚合、块编号、平局裁决都是确定性的；唯一使用伪随机数的位置是同票平局，
由 `seed` 驱动，固定种子可复现。可用 `-reread` 验证“输出可被本服务重新
读入并得到相同分相”。

## 文档

- 输入/输出文件格式：[docs/formats.md](docs/formats.md)
- 算法语义与冲突裁决策略：[docs/algorithm.md](docs/algorithm.md)
- 全部配置项：[docs/config.md](docs/config.md)
- 规模边界与资源上限：[docs/scale.md](docs/scale.md)

## 退出码

- `0`：成功（包括 `-reread` 校验通过）。
- `1`：输入解析错误、配置非法，或 `-reread` 校验失败。

非致命问题（未知位点的读段观测、非 ref/alt 等位基因、读段内部矛盾等）
不导致失败，以 `notice:` 打印到标准错误，并在统计中计数，绝不静默丢弃。
