# RPC项目统一文献列表

更新：2026-10-30。包括前期已引用和本轮实际采用的来源。

总计15项：论文10篇、官方手册2项、厂家资料2项、教学材料1项。PCB设计与用户实测输入不计入文献数量。

优先级：本装置实测/PCB/明确确认条件 > R013及其北大RPC相关文献R014/R015 > 其他RPC文献 > 扫描假设。具体参数取舍见PKU_PARAMETER_PRIORITY.md。

检索命中但未采用的来源不计入。旧代码注释中的百科链接不作为本轮电学参数的依据。



## [R001] Garfield++ User Guide

Garfield++ Collaboration，2025.1。类别：manual。

来源：[https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf](https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf)

定位：Ch. 5 Tracks, Ch. 7 Signals; Getting started weighting-field example。

用途：Heed/Magboltz接口、逐电极Ramo信号及成形函数。

使用位置：PROGRAM_STRUCTURE_AND_PHYSICS.md; scripts/electrode_signals.py; scripts/shape_strip_waveforms.py; scripts/plot_single_hit_four_ends.py; scripts/plot_shared_strip_four_ends.py; scripts/trigger_time_spectra.py; scripts/muon_trigger_time_spectra.py; TRIGGER20_READOUT.md; 第一部分报告。



## [R002] Detector physics and simulation of resistive plate chambers

W. Riegler, C. Lippmann, R. Veenhof，2003。类别：paper。

来源：[https://doi.org/10.1016/S0168-9002(03)00337-1](https://doi.org/10.1016/S0168-9002(03)00337-1)

定位：NIM A 500, 144–162。

用途：RPC电离、雪崩统计与空间电荷背景。

使用位置：PROGRAM_STRUCTURE_AND_PHYSICS.md; scripts/muon_waveform_reference.cc; 第一部分报告。



## [R003] Interfacing Geant4, Garfield++ and Degrad for the Simulation of Gaseous Detectors

D. Pfeiffer et al.，2019。类别：paper。

来源：[https://doi.org/10.1016/j.nima.2019.04.110](https://doi.org/10.1016/j.nima.2019.04.110)

定位：arXiv:1806.05880。

用途：输运/电离模块边界及避免重复计数。

使用位置：PROGRAM_STRUCTURE_AND_PHYSICS.md; scripts/muon_waveform_reference.cc; 第一部分报告。



## [R004] Muon simulations for Super-Kamiokande, KamLAND, and CHOOZ

A. Tang, G. Horton-Smith, V. A. Kudryavtsev, A. Tonazzo，2006。类别：paper。

来源：[https://doi.org/10.1103/PhysRevD.74.053007](https://doi.org/10.1103/PhysRevD.74.053007)

定位：Phys. Rev. D 74, 053007。

用途：既有Gaisser–Tang源的物理来源；非本轮新参数。

使用位置：先前源比较报告。



## [R005] Physics List Guide: Electromagnetic physics constructors

Geant4 Collaboration，online。类别：manual。

来源：[https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/electromagnetic/emphyslist.html](https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/electromagnetic/emphyslist.html)

定位：EM physics constructors。

用途：解释EM option4职责；实际程序版本以运行记录为准。

使用位置：第一部分报告。



## [R006] Currents Induced by Electron Motion

S. Ramo，1939。类别：paper。

来源：[https://doi.org/10.1109/JRPROC.1939.228757](https://doi.org/10.1109/JRPROC.1939.228757)

定位：Proc. IRE 27(9), 584–585。

用途：感应电流基本关系；符号约定同时见R001。

使用位置：scripts/electrode_signals.py。



## [R007] Series 1507 10-Tap SIP Delay Line

Data Delay Devices Inc.，1997。类别：datasheet。

来源：[https://pdf.dzsc.com/20090603/200903200838338978.pdf](https://pdf.dzsc.com/20090603/200903200838338978.pdf)

定位：Doc #97025；原厂资料转存。

用途：1507-40A每块40 ns、每抽头4 ns；非实测传输函数。

使用位置：scripts/pcb_readout_topology.py; scripts/strip_delay_readout.py; scripts/plot_strip_waveform_demo.py; scripts/shape_strip_waveforms.py; scripts/plot_single_hit_four_ends.py; scripts/plot_shared_strip_four_ends.py; scripts/trigger_time_spectra.py; scripts/muon_trigger_time_spectra.py; TRIGGER20_READOUT.md; config/rpc_response.yaml。



## [R008] Electric fields, weighting fields, signals and charge diffusion in detectors including resistive materials

W. Riegler，2016。类别：paper。

来源：[https://arxiv.org/abs/1602.07949](https://arxiv.org/abs/1602.07949)

定位：Sec. 3.2–3.3, 4.1–4.3；PDF Sec.4.1, Eq.(4.10)。

用途：多介质/阻性材料权重场；典型玻璃εr=8、体电阻率量级10^10–10^12 Ω cm。

使用位置：STRIP_READOUT_STATUS.md; material scan。



## [R009] Signal coupling and signal integrity in multi-strip Resistive Plate Chambers used for timing applications

D. Gonzalez-Diaz, H. Chen, Y. Wang，2011。类别：paper。

来源：[https://doi.org/10.1016/j.nima.2011.05.039](https://doi.org/10.1016/j.nima.2011.05.039)

定位：NIM A 648, 52–72；arXiv:1102.1389。

用途：模态色散、耦合和传输完整性；说明纯延迟模型适用边界，未挪用其传输参数。

使用位置：STRIP_READOUT_STATUS.md。



## [R010] Studies on Electrical Properties of Resistive Plate Chamber (RPC)

S. Das, J. Datta, N. Majumdar, S. Mukhopadhyay，2022。类别：paper。

来源：[https://arxiv.org/abs/2204.12986](https://arxiv.org/abs/2204.12986)

定位：v2, Sec.3.1 Table 3；JINST 17 P09041。

用途：该文模型：玻璃εr=10、ρ=10^12 Ω cm；石墨Rs=300 kΩ/□、εr=10、厚0.1 mm；非本RPC实测。

使用位置：config/strip_material_scan.json; material scan。



## [R011] PET & PEN Films for Capacitors

DuPont Teijin Films / Mylar Specialty Films，2024 brochure。类别：datasheet。

来源：[https://www.dupontteijinfilms.com/wp-content/uploads/2024/02/Mylar-Specialty-Films-Product-Capacitors-Brochure.pdf](https://www.dupontteijinfilms.com/wp-content/uploads/2024/02/Mylar-Specialty-Films-Product-Capacitors-Brochure.pdf)

定位：Typical properties table, 25°C and 1 kHz。

用途：Mylar PET εr=3.25；用作PET代理，非本膜材料测量或全频介电谱。

使用位置：config/strip_material_scan.json; material scan。



## [R012] Design and construction of RPC detectors

G. Aielli，DRD1 school slides。类别：teaching。

来源：[https://indico.cern.ch/event/1384298/contributions/6070891/attachments/2978847/5244587/RPC%20construction_DRD1%20school.pdf](https://indico.cern.ch/event/1384298/contributions/6070891/attachments/2978847/5244587/RPC%20construction_DRD1%20school.pdf)

定位：slide 21 Resistive Plate – ohmic contact。

用途：石墨方阻约100 kΩ/□及快速信号透明性的时间尺度解释；教学量级。

使用位置：config/strip_material_scan.json; material scan。



## [R013] A proposed PKU-Muon experiment for muon tomography and dark matter search

Xudong Yu et al.，2024。类别：paper。

来源：[https://arxiv.org/abs/2402.13483](https://arxiv.org/abs/2402.13483)

定位：用户本地PDF为v2, 2024-03-02；PDF p.4 Fig.6及正文，p.14 refs.[63–66]；期刊版DOI 10.1103/PhysRevD.110.016017（不能混同版本页码）。

用途：优先装置谱系：float glass、2 mm气隙、100 μm PET、4 ns延迟；原型2.6 mm玻璃与当前2.7 mm不同；文中未来混合气不覆盖当前纯R134a。

使用位置：PKU_PARAMETER_PRIORITY.md; config/strip_material_scan.json; STRIP_CONNECTION_LOGIC.md; scripts/strip_delay_readout.py; scripts/plot_strip_waveform_demo.py; scripts/shape_strip_waveforms.py; scripts/plot_single_hit_four_ends.py; scripts/plot_shared_strip_four_ends.py; scripts/trigger_time_spectra.py; scripts/muon_trigger_time_spectra.py; TRIGGER20_READOUT.md。



## [R014] Study of spatial resolution properties of a glass RPC

Qite Li, Yanlin Ye, Chao Wen, Wei Ji, Yushou Song, Rongrong Ma, Chen Zhou, Yucheng Ge, Hongtao Liu，2012。类别：paper。

来源：[https://doi.org/10.1016/j.nima.2011.10.006](https://doi.org/10.1016/j.nima.2011.10.006)

定位：NIM A 663,22–25；Sec.2,p.23；作者单位公开全文 https://researchgroups.pku.edu.cn/_tsf/00/0E/mYfayyMfyQnq.pdf；R013引用[63]。

用途：石墨方阻约10^7 Ω/□；同系列原型与前端历史信息。300 MHz与放大器型号不唯一确定传递函数。

使用位置：PKU_PARAMETER_PRIORITY.md; config/strip_material_scan.json; STRIP_CONNECTION_LOGIC.md; scripts/strip_delay_readout.py; scripts/plot_strip_waveform_demo.py; scripts/shape_strip_waveforms.py; scripts/plot_single_hit_four_ends.py; scripts/plot_shared_strip_four_ends.py; scripts/trigger_time_spectra.py; scripts/muon_trigger_time_spectra.py; TRIGGER20_READOUT.md。



## [R015] A sub-millimeter spatial resolution achieved by a large sized glass RPC

Qite Li, Yanlin Ye, Wei Ji, Chao Wen, Hongtao Liu, Yucheng Ge，2013。类别：paper。

来源：[https://doi.org/10.1088/1674-1137/37/1/016002](https://doi.org/10.1088/1674-1137/37/1/016002)

定位：Chinese Physics C 37,016002；Sec.2,p.016002-2；R013引用[64]。

用途：再次核对10^7 Ω/□、100 μm PET、1507-40A及4 ns抽头；约4mm电荷轮廓FWHM用作高斯共享演示的宽度参考，非本PCB拟合；所报坐标FWHM与本项目σ分开。

使用位置：PKU_PARAMETER_PRIORITY.md; config/strip_material_scan.json; STRIP_CONNECTION_LOGIC.md; scripts/strip_delay_readout.py; scripts/plot_strip_waveform_demo.py; scripts/shape_strip_waveforms.py; scripts/plot_single_hit_four_ends.py; scripts/plot_shared_strip_four_ends.py; scripts/trigger_time_spectra.py; scripts/muon_trigger_time_spectra.py; TRIGGER20_READOUT.md。



## 非文献输入

PCB：用户提供文件，原路径及SHA256见 config/pcb_geometry.json。FR4设计εr=4.8，不作为文献测量值。

几何厚度、2 mm气隙及低计数率稳态电压近似来自原配置和用户确认。0.749 mm为空间高斯σ实测输入。

扫描点是研究预设，不代表材料不确定度或文献置信区间。