# RPC 响应参考库与坐标基线

模块职责、数据流和物理公式见[程序结构与物理原理说明](PROGRAM_STRUCTURE_AND_PHYSICS.md)。

当前阶段由用户确认：先交付线性雪崩参考库和坐标基线；前端阈值、效率、电荷谱未提供，不能据非零雪崩预测实际探测成功。

## 已确认条件

| 参数 | 当前值及来源 |
|---|---|
| 气体 | config/newrpc_material.yaml 的纯 C2H2F4（R-134a） |
| 温度、压力 | 用户指定 298.15 K、760 Torr |
| 外电极电压 | DetectorConstruction.cc 的 10,100 V |
| 有效气隙 | newrpc.yaml 的 2 mm，不包括两侧 2.7 mm 玻璃 |
| 场近似 | 用户确认低计数率直流稳态，忽略玻璃压降，气隙场 50,500 V/cm |
| 孤立命中核心分辨率 | X/Y 均取 0.749 mm 高斯 sigma；默认作为总误差，不再加雪崩扩散展宽 |
| 延迟线 | 用户确认每方向 8 个 1507-40A 完整串联 |
| 名义延迟 | 单块40 ns、10个4 ns间隔；整链320 ns，80间隔对应81个读出位置 |

config/rpc_response.yaml 是物理假设与待标定参数清单，当前命令行程序仍使用显式参数；该 YAML 不是自动驱动全部程序的运行配置。未知值保留 null。

## 程序与数据协议

* rpc_make_gas：Magboltz 气体输运表；指定场强扫描范围和统计量。
* rpc_make_response：从 .gas 取指定工作场的输运参数，导出单位明确的 .rpc 和 ROOT 文件。范围外查询拒绝外推。
* rpc_make_kernel：按起始漂移距离生成单电子线性出生/死亡雪崩样本，保存理想平板电子感应电荷及有效性。
* rpc_query_kernel / include/RpcKernel.hh：一次载入 ROOT 样本，内存中按最近深度格抽样；不依赖运行时 Garfield++。
* rpc_digitize：读取 Geant4 rpcSegments，生成经验坐标基线及名义四端时间，不调用 Garfield++。

所有新运行时数据使用 **mm、ns、eV**。气体文件内部遵循 Garfield++ 的 cm、ns、V/cm；转换在导出器完成。新参考 ROOT 格式与旧 rpc_cluster_params.root 的 clusterParams 半径表不同，不可互换。

### Geant4 新输出

原有 tree/Edeps/Tracks 保留。新增平坦 TTree `rpcSegments`，在旧四层触发筛选之前记录真实2 mm气隙内的步，包括零沉积步：

`event, layer, track, parent, pid, x0_mm, y0_mm, z0_mm, x1_mm, y1_mm, z1_mm, t0_ns, t1_ns, edep_ev, nonion_ev, kinetic_ev, length_mm, gap_center_z_mm, gap_mm`。

由现有水平分层结构的读出板中心及 film/electrode/glass/gas 厚度计算气隙中心。当前仅支持现有水平对齐布局；若修改为倾斜、镜像或不同堆叠，必须先扩展为几何实例局部坐标变换。每个文件的 event 是生成事例编号；不同文件合并时须加文件来源键，不能直接按 event 把事例混合。

新记录保留次级粒子来源、步端点和时间，避免旧 Edeps 按种类/过程合并后无法恢复轨迹的问题。Geant4 多重散射仍受其步内凝聚近似约束。新导出不修订旧 Tracks 分支里历史的 E 字段混用逻辑。

### 坐标基线

在每个事例、每层内，以用户显式给出的时间窗组成非滑动时间组；按局部 `edep_ev - nonion_ev` 权重计算步中点重心。它是基准重心近似，不是条带波形重建。

读出结果保存为 `rpcHits`：event/layer、x_mm/y_mm、t_ns（组内最早沉积中点时间）、ionizing_edep_ev、n_segments/n_tracks、ambiguous、valid_coordinate，以及 nominal_txl_ns/txr_ns/tyd_ns/tyu_ns。

默认多轨迹组（包括次级轨迹）标 ambiguous=true、valid_coordinate=false，坐标/端点时间写 NaN；不把 sigma 当双命中合并距离。`--allow-multi-track` 仅显式开启能量重心近似，仍保留 ambiguous=true，不能当作配对已解决。单轨迹长程弯折也只按重心近似处理。

当前每方向的顺序抽头映射采用名义假设：由现有读出几何得到 X 周期3.457 mm、Y 周期3.463 mm，相邻位置时间差变化8 ns，所以 ax=3.457/8、ay=3.463/8 mm/ns。X读出条中心相对板中心偏移1.062 mm，Y偏移0。外侧条中心跨度分别276.56和277.04 mm；不能拿板长300 mm直接除以总延迟。通道顺序和方向最终应以位置扫描标定替换。

为满足实测空间 sigma，构造等方差、相互独立的端点抖动：sigma_end = sigma_position / (sqrt(2)*abs(a))；再由时间差恢复坐标。这里只是符合空间宽度的一种实现，不唯一。名义四端时间不用于预测绝对时间分辨率、电子学共模抖动、阈值时移或双脉冲能力。

没有延迟线波形、死时间与配对标定时，不进行跨 Geant4 event 的真实到达流混合；320 ns传播延迟也不能直接当作死时间或双脉冲分辨时间。

### 单电子参考库

使用常系数 Townsend 电离与附着的线性出生/死亡过程。每一小段抽样其精确人口转移分布，电子感应电荷采用理想平板 Ramo 梯形积分。包含附着/消失和增益涨落；不模拟离子慢尾、真实条带权重场、空间电荷或电子学阈值。

MAX_ELECTRONS 是数值/线性模型有效性保护，不是已标定的饱和电荷；越界样本 valid=false。查询不得丢弃该样本后重抽，否则会人为压低大雪崩比例。最近格查询的起始深度误差至多半格，需按目标精度检查网格收敛。引用的 .rpc 参数快照内嵌于 ROOT，便于追溯。

该参考库暂不用于改变坐标基线的重心权重。将 Geant4 局部沉积转换为微观电子数时仍需 W、簇统计和 Geant4/Heed 职责验证；不能对 Geant4 已输运的 δ 电子再重复加完整 Heed 级联。

## Ubuntu 构建和运行

共享路径为 `/home/randomuon/muon/PKMUON_2024-HCR`。安装路径可用 GARFIELD_INSTALL、GEANT4_INSTALL 覆盖。

```bash
cd /home/randomuon/muon/PKMUON_2024-HCR
source rpc_response/scripts/env_ubuntu.sh
bash rpc_response/scripts/build_ubuntu.sh
mkdir -p rpc_response/results
```

气体表为离线步骤，生产事件不自动生成；ncoll 的单位是 Magboltz 的 10^7 次碰撞。首轮扫描 ncoll=1 为探索性统计，不代表精密输运标定。

```bash
rpc_response/build-ubuntu/rpc_make_gas rpc_response/results/r134a_298K_760Torr_scan.gas 10000 60000 6 1
rpc_response/build-ubuntu/rpc_make_response rpc_response/results/r134a_298K_760Torr_scan.gas 2 10100 rpc_response/results/rpc_transport_v1.rpc rpc_response/results/rpc_transport_v1.root
rpc_response/build-ubuntu/rpc_make_kernel rpc_response/results/rpc_transport_v1.rpc rpc_response/results/rpc_linear_reference_v1.root 40 1000 20261018 10000000 0.01 --linear-reference
rpc_response/build-ubuntu/rpc_query_kernel rpc_response/results/rpc_linear_reference_v1.root 1.0 1000000 20261018
```

响应导出允许在场强网格内部插值；当前50.5 kV/cm在扫描表内，不是单独的精密工作点计算。

运行固定种子的200事件气隙检查（MUPOS_RPC_PENCIL_BEAM只改变验证时的初级源，未设置时保持原Gaisser-Tang源）：

```bash
cd build-rpc
MUPOS_RPC_PENCIL_BEAM=1 ./muPos ../rpc_response/tests/smoke.mac
cd ..
rpc_response/build-ubuntu/rpc_root_check segments rpc_response/results/geant4_smoke.root
rpc_response/build-ubuntu/rpc_digitize rpc_response/results/geant4_smoke.root rpc_response/results/geant4_baseline.root 10 20261018
```

上面的10 ns只用于演示分组流程，不是测得的双脉冲分辨能力。新工具保护已存在输出，数字化输出先写 `.partial` 再改名；报错后的 partial 是未完成文件，不能作为正式结果。Geant4 的原有输出行为仍为覆盖写，重跑前应另设输出路径。

## 验证与遗留路径

内核验证覆盖纯出生/纯附着/临界/净附着的均值、Ramo积分、延迟线无噪声往返、时间分组、0.749 mm高斯宽度及越界保护。ROOT测试覆盖已知重心的20000事例；Geant4测试检查200事例、四层、所有步均在2 mm气隙、存在次级轨迹。

旧 analysis.C 的 SIGMA_CM 数值误用于 Geant4 毫米坐标，现已改为0.749 mm；旧 LUT 半径按其声明的cm转换为输出mm。旧宏的2sigma合并和选最大能量簇属于遗留算法，新流程不使用。旧程序运行时自动启动 Garfield 的逻辑改为仅显式设置 RPC_LEGACY_AUTOGENERATE 才启用，避免多进程重复计算。

Geant4材料配置仍保留原指定密度0.00425 g/cm3；它不是按298.15 K重新计算的理想气体密度。未来做定量电离/效率比较时，需要将其与Magboltz数密度一致化，或确认该密度为实测值。当前参考库与经验坐标基线不构成已标定的完整RPC效率/电荷模型。

## 方法来源

* [Garfield++ 手册](https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf)：气体输运系数、单位与Shockley–Ramo感应信号。
* [Riegler–Lippmann–Veenhof](https://doi.org/10.1016/S0168-9002(03)00337-1)：初级电离、统计雪崩与电子学的分层思路。当前代码是受其结构启发的基础线性模型，未声称复现其完整RPC模型。
* [Geant4/Heed接口研究](https://arxiv.org/abs/1806.05880)：电离与次级粒子职责划分。
* [1507原厂数据表存档](https://pdf.dzsc.com/20090603/200903200838338978.pdf)：1507-40A名义总延迟40ns、每抽头4ns；制造公差不当作逐事件随机时间抖动。

没有复制第三方模拟源码；后续若直接复用代码，需要单独核对许可证。

## 缪子标准样本与分布图

新增 rpc_make_muon_standard，以 **1 GeV、垂直入射的mu-** 建立当前标准工作点。不是对全部能量与入射角通用的响应表。
独立Heed拥有本样本的完整初级电离及δ电子级联，不叠加Geant4电离。样本位置在X/Y各[-100,100] mm均匀抽取，气隙2 mm；使用Heed电子的实际出生位置/时间查询单电子库。

```bash
# 在Ubuntu项目根目录，先加载env_ubuntu.sh
rpc_response/build-ubuntu/rpc_make_muon_standard \
  rpc_response/results/rpc_transport_v1.rpc \
  rpc_response/results/rpc_linear_reference_v1.root \
  rpc_response/results/r134a_298K_760Torr_scan.gas \
  rpc_response/results/muon_1GeV_normal_v1 10000 1000000000 20261019
```

输出muonStandard ROOT树、events/electrons/clusters/kernel/transport CSV和运行元数据JSON。现有输出受保护，重跑需换前缀。
任一电子查询无效时，整个缪子的总线性电荷设NaN，不能用其余有效电子的电荷和冒充完整响应。
电子出生深度对应向-z收集面的漂移距离；creation+distance/velocity仅为忽略扩散的名义到达时间，不能把它当作实际阈值时间。

坐标仍用真实穿越位置与实测总sigma构造名义四端时间；不是对Heed电离云中心再增加一份sigma。这样避免将已包含于实测sigma的展宽重复计入。
绝对端点时间不包含已标定的雪崩起始抖动或电子学共模时间误差。所有坐标残差是在假定存在孤立读出信号的条件下生成，不预测触发效率。

在安装numpy/matplotlib的环境执行：

```bash
python rpc_response/scripts/plot_muon_standard.py rpc_response/results/muon_1GeV_normal_v1
```

生成四组PNG、单页PDF、合并PDF和summary.json：电离统计、线性参考库有效范围、名义延迟线与坐标残差、气体输运扫描。
电荷图明确筛选valid样本；输运扫描的连线仅辅助观察，50.5 kV/cm为表内插值；图中无效比例不是探测器效率。

## 多能量标准库与50 ns事例比较（2026-10-20）

**后续更新：** 已增加用户确认的空间多位置近似，并完成标准库插值交叉验证。最新入口为 `results/comparison_20261020/spatial_interpolation/spatial_interpolation_report.html`；该报告替代下文“多轨迹一律拒绝”作为唯一策略的描述，旧结果保留对照。

```bash
# Ubuntu：复用已生成的10万次平面源结果；默认5σ，扫描另含3σ与10σ
bash rpc_response/scripts/run_spatial_interpolation.sh
# Python（需numpy/pandas/scipy/matplotlib）：全量守恒检查、留一及独立中间点验证与图表
python rpc_response/scripts/analyze_spatial_interpolation.py rpc_response/results/comparison_20261020
# 独立使用插值查询；范围0.1—100 GeV，输出整条事例抽样
python rpc_response/scripts/response_interpolation.py \
  rpc_response/results/comparison_20261020/energies 2 /tmp/interpolated_2GeV.csv
```

新ROOT结果每个空间位置一行，可按 `(event,layer,window_index)`取同组所有坐标。
`n_positions`为同组位置数，`n_tracks`为该位置内轨迹数，`window_n_tracks`为同组轨迹数；`centroid_x_mm/y_mm`保留噪声前中心。
`valid_coordinate`只表示在本模型假设下有效；`ideal_spatial_model=true`、`readout_pairing_validated=false`明确说明空间近似与实际延迟线配对的区别。
小于或等于所选分离门限的中心连通归并，远距离分量分别展宽；一条延展轨迹仍按一个中心处理，幅度阈值未标定。

插值默认在log(E)上对两个相邻经验分布取凸组合，通过抽取完整事例保持整数计数及内部相关性；没有混合两条事例的整数电子数，也不静默删除无效电荷。
低能电离谷附近的留一误差较大，不能移除0.3 GeV标准点；垂直入射库不能凭能量插值推广到任意角度或粒子种类。

本轮结果入口为 `results/comparison_20261020/comparison_report.html`，中文说明另存 `COMPARISON_REPORT.md`，全部科学图表合并为 `figures/all_comparison_figures.pdf`。
已生成0.1、0.3、1、3、10、100 GeV六个标准点，每点10000个垂直缪子。复现入口为 `scripts/run_energy_scan.sh`。

Geant4运行保留原−Y球面源，并新增 `MUPOS_GAISSER_PLANE=1` 平面模式：负Z端RPC外侧1 mm、400×400 mm平面、朝+Z，联合拒绝采样要求直线穿过四层真实气隙。接受后的源谱具有几何条件，不能直接当作绝对宇宙线通量。`rpcPrimaries`记录实际源坐标、方向、模式和提议次数；`rpcSegments.root_track`用于区分次级轨迹与独立初级入射。运行入口为 `scripts/run_gaisser_comparison.sh plane`。

50 ns已确认为事例分组/符合窗；名义320 ns延迟线传播后仍允许采集四端信号，不能将50 ns作为原始四端时间戳门宽。实际双脉冲分辨时间仍未标定。
`scripts/analyze_comparison.py`生成能量、源几何、坐标残差、自然次级轨迹、受控独立双缪子叠加与理想端点配对图；`scripts/build_comparison_report.py`生成带图中文报告。

当前默认数字化仍对多轨迹组输出歧义，不输出多个坐标。理想配对原型独立于默认数字化，假设所有端点脉冲均已提取，不能当作电子学已实现的双命中效率。
线性雪崩参考在当前工作点有大量总电荷无效事件，因此本轮坐标比较使用真实气隙沉积与实测sigma基线，没有声称完成Garfield++电荷驱动的坐标重建。报告同时保留覆盖率与长尾，避免只比较有效单轨迹残差。
