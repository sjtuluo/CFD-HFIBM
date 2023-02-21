CFD Model HydroFlow Immersed Boundary Method
===
本项目是计算流体力学程序HydroFlow-IBM的数值模型框架，致力于完成课题项目**波流-结构物耦合的浸没边界法数值模拟研究**。  
程序主框架HydroFlow由上海交通大学船舶海洋与建筑工程学院张景新副教授课题组开发。HydroFlow数值模型基于非结构网格、二阶TVD格式的有限体积法建立，主要应用于具有自由表面流动的物理问题数值模拟，如近海岸洋流、河流、湖泊等自然地表水系与环境相互影响的水动力学特征研究。  
该数值模型使用垂向坐标变换法模拟流体自由表面位置，并使用半隐式格式进行数值离散。计算程序配有雷诺平均方程模型（Reynolds-averaged Navier-Stokes，RANS）、大涡模拟（Large Eddy Simulation，LES）、分离涡模型（Detached Eddy Simulations，DES）三种湍流计算模式可供选择，同时开发离散元多孔介质模型（Discrete Element Method，DEM）、浸没边界法模型（Immersed Boundary Method，IBM）、两相流模型（Two Phase Flow Model）等功能模块，可以应用于溃坝模拟、泥沙输运模拟、水生植物群落运动模拟、波流结构物耦合模拟等项目研究。  
本项目在HydroFlow数值模型基础上进一步完成浸没边界法模型（IBM）的开发与应用。  
更多数值模型详细内容请参考论文[An efficient 3D non-hydrostatic model for simulating near-shore breaking waves](https://github.com/sjtuluo/CFD-HFIBM/blob/master/Documents/HydroFlow%20Model.pdf)。  

## 程序预设文件
程序运行之前，需读入计算网格几何信息文件并修改边界条件，同时应设置计算模型的相关参数。预设文件简要信息介绍请参考以下表格。
| 文件名 | 文件内容 |
| ------ | ----------- | 
| OCERM_INF | 计算程序头文件，计算信息与全局变量设置 |
| OCERM.GRD | 计算网格文件，包含垂向分层信息、网格节点位置与节点编号等内容 |
| OCERM.CUV | 计算网格文件，包括网格单元拓扑关系、网格边与局部坐标信息等内容 |
| infl.QBC | 入口流量边界条件文件，按照入口单元水平位置与垂向分层设置流量 |
| outl.EBC | 出口水位边界条件文件，按照出口单元位置设置水位 |
| VIS.QBC | 入口湍流边界条件文件，按照入口单元水平位置与垂向分层设置相应湍流参数 |
| VIS.EBC | 出口湍流边界条件文件，按照出口单元水平位置与垂向分层设置相应湍流参数 |
| Gauge_XY.DAT | 监测点位置文件，按照监测点个数与对应位置逐一设置 |

  
## 计算模块简介
程序包含基础模块较多，这里仅对一个完整时间步内运行的基本计算模块作简要介绍。完整的程序模块说明请参考[模型结构及变量说明](https://github.com/sjtuluo/CFD-HFIBM/blob/master/Documents/%E6%A8%A1%E5%9E%8B%E7%BB%93%E6%9E%84%E5%8F%8A%E5%8F%98%E9%87%8F%E8%AF%B4%E6%98%8E.pdf)，关于模型算法的介绍与程序实现请参考[模型算法简介](https://github.com/sjtuluo/CFD-HFIBM/blob/master/Documents/%E6%A8%A1%E5%9E%8B%E7%AE%97%E6%B3%95%E7%AE%80%E4%BB%8B.pdf)与[程序模块注解](https://github.com/sjtuluo/CFD-HFIBM/tree/master/Documents/%E7%A8%8B%E5%BA%8F%E6%A8%A1%E5%9D%97%E6%B3%A8%E8%A7%A3)。  
基本计算模块简明介绍请参考以下表格，模块说明顺序按照程序运行顺序给出。  
- 计算信息读入部分  

| 计算模块名 | 计算模块功能简介 |
| ------ | ----------- | 
| SETDOM | 计算网格信息读取，初始化计算域相关几何信息 |
| BCDATA | 边界条件信息读取，初始化相应边界变量 |
| ZEROES | 全局变量初始化 |

- 时间递进计算部分   

| 计算模块名 | 计算模块功能简介 |
| ------ | ----------- | 
| BCOND | 边界条件计算模块，通过程序参数控制不同边界条件模式 |
| SUBGRIDH / SUBGRIDV | RANS湍流模型模块（预设），模块内可使用S-A湍流模型与SST湍流模型计算涡粘系数，使用不同湍流模型需对应不同边界条件设置 |
| ADVU / ADVV / ADVW | 动量方程离散与求解模块，对流扩散项、源项与显式水位梯度项计算 |
| TVDSCHEMEH / TVDSCHEMEV| TVD（Total variation diminishing）模型计算模块，动量方程对流项采用二阶TVD格式计算面通量 |
| ELTION | 水位求解模块，计算模型采用时间半隐式格式离散，通过连续性方程求解水位，水位求解方程系数矩阵为对角矩阵，应用双共轭梯度法迭代求解 |
| PROFV / PROFW | 流速求解模块，计算模型采用分步投影法求解流速，流速求解方程系数矩阵为三对角矩阵，应用追赶法迭代求解 |
| VERTVEL | 静压模型垂向速度求解模块，计算坐标变换后的垂向速度 |
| DYN | 动压模型计算模块，动压求解方程系数矩阵为对角矩阵，应用双共轭梯度法迭代求解 |
| WREAL | 物理垂向速度求解模块，计算不受坐标变换影响的物理垂向速度 |
| IBM | 浸没边界法模块，本项目研究开发方向 |
| REUV | 流场变量更新模块，完成一个时间步计算更新相应变量并计算真实物理流速用于后处理 |
| UVFN | 流场变量插值模块，计算网格边与顶点的流场变量并更新网格局部信息 |
| ARCHIVE | 数据存储模块，输出指定时间步计算结果文件，依照Tecplot数据格式生成 |

- 浸没边界法部分  

| 计算模块名 | 计算模块功能简介 |
| ------ | ----------- | 
| IBMINIT | 浸没边界法模块初始化与边界数据读取 |
| IBMALIDC / IBMALIDP | 虚拟边界识别算法，划分固体域与流体域，识别Ghost Cells |
| IBMALDF | Delta Function方法计算虚拟力源项 |
| IBMALGC | Ghost Cell方法计算虚拟力源项 |
| IBMRGD | 运动边界与刚体动力学模型 |
| IBMCLBD / IBMCLBP | 刚体碰撞检测模型 |
| IBMVISUAL| 虚拟边界识别结果储存与可视化 |
| IBMARCHIVE | 浸没边界法模块计算结果存储 |

以上为时间递进过程中单一时间步内基本计算模块简明介绍，还可根据需求应用其他模块。时间递进至设定值后计算完成，可查看相应结果文件并进行数据处理。
