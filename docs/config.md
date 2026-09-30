# 配置项

配置为可选 JSON 文件（`-config`）。缺省文件或缺省字段使用默认值；非法
取值在启动时被拒绝并返回非零退出码。完整配置会随分相结果写入 `#config`
段，因此结果文件自描述且可复现。

| 键 | 类型 | 默认值 | 合法范围 | 语义 |
|---|---|---|---|---|
| `min_edge_support` | int | `2` | `>=0` | 边连通所需最少总票数（读段票与家系票合计）；`0` 关闭深度要求 |
| `min_edge_strength` | float | `0.8` | `[0.5,1]` | 多数票最低占比；低于则断块。低于 0.5 无裁决意义故拒绝 |
| `min_block_confidence` | float | `0.95` | `(0,1]` | 块置信度低于该值标 `LOW_CONF`（不影响连/断边） |
| `block_max_error_rate` | float | `0.25` | `[0,1]` | 边预期错误率上限；超过则断块。默认 0.25 允许 2 票全一致（平滑后 0.167）连通 |
| `min_site_depth` | int | `3` | `>=0` | 有效读段覆盖数下限；低于标记 `low_cov=1`，不影响分相成败 |
| `seed` | int64 | `1` | 任意 | 仅用于等票平局的确定性裁决；固定种子可复现 |
| `enable_trio` | bool | `true` | bool | 是否启用家系孟德尔约束与父源定向 |
| `enable_reads` | bool | `true` | bool | 是否启用读段连边证据 |

示例见 `examples/config.json`。

## 阈值如何影响分块

- 连通条件三者同时满足：`min_edge_support`、`min_edge_strength`、
  `block_max_error_rate`。任一不满足即在该相邻杂合边处分块。
- `min_block_confidence` 只决定块上的 `LOW_CONF` 标记，不会拆分块。
- 关闭读段（`enable_reads=false`）且无可用家系边时，每个杂合位点都是
  `SINGLETON` 块；这是确定性的退化结果，不报错。
