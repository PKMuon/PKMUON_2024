# 条带波形20%判时与位置读出（2026-10-31）

新增1GeV源路径：`scripts/muon_waveform_reference.cc`读取既有Heed逐电子CSV，重新抽样各电子家系雪崩并累积事例电流；`scripts/muon_trigger_time_spectra.py`再执行共享、延迟、成形和20%判时。输出见 `../deliverables/Muon1GeV_Trigger20_20261031/`。10000个源事例中9974个超出线性参考界限，全量图仅为未饱和形式参考，26个界限内子样本另图，不能用任一图宣称已预测完整RPC响应。原下面5000样本均匀位置研究仍保留作为算法对照，不与1GeV实际电离事例混同。

旧经验重心入口保留对照。新增的 `scripts/digitize_strip_waveforms.py` 从四端成形波形重新提取时间并重建位置，不直接使用真值位置或端点峰时。

每端先取本端正峰值A，阈值为用户指定0.2A，取峰前上升沿首次越阈处，两采样点之间线性插值。无正信号、采集窗截断、多个显著脉冲等输出状态标志，歧义事件不输出坐标。多峰辅助判据为峰高至少0.2A且prominence至少0.05A，后者为算法诊断预设，不是实测双脉冲能力。该离线20%方法不等同于硬件CFD电路，不预测触发效率。

对每方向，`tap=(T_net1-T_net81+320 ns)/(2*4 ns)`。将tap在实际PCB通道中心间插值得到位置，再减去中心通道坐标。公共起始时间抵消。多条带波形形状引入的判时偏差未自动校正。

运行：

```
python rpc_response/scripts/digitize_strip_waveforms.py four_end_waveforms.csv positions.csv
python rpc_response/scripts/trigger_time_spectra.py --events 5000
python rpc_response/tests/test_trigger_readout.py
```

依赖：numpy、pandas、scipy、matplotlib。样例雪崩CSV和PCB网络/铜几何已放入examples、config，可随仓库获取，无需.part1_build目录。完整Heed/Garfield参考库仍需按既有工具生成。

本轮5000事例是均匀独立x,y∈[-100,100] mm的模型样本，不是新运行的Geant4或Gaisser–Tang宇宙线。采用同一已存在单电子线性雪崩模板、4 mm FWHM高斯几何共享、4 ns抽头、5 ns单极成形、t0=0，不加噪声、实测σ或独立雪崩涨落。

二维谱横轴为该方向网络1端20%时间，纵轴为网络81端20%时间，单位ns，颜色为每2 ns×2 ns箱的事例数。固定t0导致反相关窄带。时间和均值约335.56 ns，多出的约15.56 ns包含两端判时/模板时移，不能把峰前20%时间直接当传播延迟或认定总链延迟改变。

该无噪声模型下X/Y位置RMS误差约0.2735/0.2736 mm，最大绝对偏差约0.392 mm。它是规则条带共享和20%判时造成的模型偏差，不是RPC空间分辨率预测。

对64个相同事件将采样从0.1 ns细化到0.02 ns，四端时间最大变化约0.00140 ns。判时基本测试及原条带积分守恒测试通过。

结果在 `../deliverables/Trigger20_Spectra_20261031/`。包含四端时间与坐标CSV、二维谱PNG/PDF、位置偏差图和单命中20%标记示例。

引用：R001成形函数；R007/R014/R015延迟链和历史电荷轮廓。20%阈值为本次用户指定，无新外部文献。完整总表为REFERENCES.md。
