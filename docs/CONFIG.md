# 配置项与规模边界

配置通过 `-config` 指向的 JSON 文件提供。所有字段均可省略，省略或为零值时使用
下表默认值。实际生效的完整配置会原样写回结果文件的 `config` 字段。

| 字段 | 类型 | 默认值 | 合法范围 | 含义 |
| --- | --- | --- | --- | --- |
| `min_edge_support` | int | `1` | `>= 1` | 连接两个相邻杂合位点所需的独立读段投票下限 |
| `min_reads_per_site` | int | `1` | `>= 1` | 位点读段覆盖下限，低于则标记 `low_coverage` |
| `min_block_confidence` | float | `0.90` | `[0.5, 1.0]` | 接受读段边所需的后验置信度下限 |
| `read_error_rate` | float | `0.01` | `(0, 0.5)` | 单观测错误概率，用于后验置信度与裁决 |
| `seed_enabled` | bool | `false` | — | 是否用随机种子裁决同票边 |
| `random_seed` | int | `0` | 任意 | 平局裁决随机种子（仅 `seed_enabled=true` 时生效） |

说明：

- `min_edge_support` 与 `min_block_confidence` 是“断开分相块”的唯一阈值；二者都满足
  才连接，否则在该相邻位点处断块，行为确定。
- 当 `seed_enabled=false`（默认）时，同票边按参与读段中最小的 `read_id` 投票裁决，
  结果完全确定且与随机数无关。
- 非法配置（范围越界等）会在运行前报错退出，不会回退到隐式取值。

示例：

```json
{
  "min_edge_support": 2,
  "min_reads_per_site": 3,
  "min_block_confidence": 0.95,
  "read_error_rate": 0.01,
  "seed_enabled": true,
  "random_seed": 42
}
```

## 声明的处理规模边界

在 Apple Silicon / arm64 与 Go 1.25+（amd64 等价）上，本服务声明支持以下规模：

| 维度 | 边界 |
| --- | --- |
| 变异位点数 | 10,000 |
| 样本数 | 100 |
| 读段观测总数 | 2,000,000 |

在该边界规模下，由自动化规模测试 `TestScaleBoundary` 守护，资源上限为：

| 指标 | 上限 | 实测参考 |
| --- | --- | --- |
| 端到端运行时间（解析 + 分相 + 写结果） | 30 秒 | 约 2.5–4 秒 |
| 运行后保留堆内存（GC 后） | 512 MiB | 约 256 MiB |

规模测试可用以下命令运行（约数秒）：

```sh
go test -run TestScaleBoundary -v .
```

短期 `go test` 冒烟可跳过规模测试：

```sh
go test -short ./...
```

资源随规模平稳增长：解析为流式扫描，内存不随文件行数线性堆积；每样本证据以稀疏边
（至多 `杂合位点数-1` 条）与近线性并查集表示，不随读段数量保留逐观测中间结构。
