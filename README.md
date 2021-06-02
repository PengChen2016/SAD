# SAD  
**the Simple Analysis of the Driver in RF hydrogen ion source**  



## get collision cross section
https://nl.lxcat.net/data/set_type.php -> type: Scattering crossing sections -> database: all -> specA: Ground states: e -> specB: Ground states: Ar -> groups: Elastic, Ionization -> processes: Biagi-v7.1, Morgan -> compare and download -> put each data set into a single file, record the reference at the first line -> modify the get_Xsec.m  
Cross sections of e-H2, e-Ar has been included in the ./packages/physic_tools.  

## compare data
output data by get_output_json(), then compare file by arbitrary text comparator, such as plug-in of code editor

## History before version control by Git
v190111 等效媒质模型与变压器模型初版自左晨，赵鹏  
190916 by陈鹏 进行等效媒质模型代码整理  
v200413 by陈鹏 完成等效媒质模型与变压器模型检查，单点计算可靠  
以2018Jain中ELISE case做benchmark：见main_v200413.m和result_of_main_v200413.md  
CHARLIE case：见main_v200413_2.m和result_of_main_v200413_2.md  
1.随机加热频率vst根据1995Vahedia修改  
2.变压器模型中二次线圈理论计算式进行了长冈系数校正；互感表达式考虑了线圈电感。校正后在ELISE case中与2018Jainb结果相近，在CHARLIE中比2018Jainb结果更优  
  
v200418 by陈鹏 进行代码修改  
有损介质中波长公式，考虑电导率影响;考虑有限半径几何效应的等效集肤深度;M考虑Lp  
分析见result_of_main_v200418.md，并提供新的ELISE case benchmark  
  
v181211 解析电模型初版自李增山。  
v200420 by陈鹏 将解析电模型加入main。  
与LZS结果做对比，有较大区别，结果见result_of_main_v200420_1.md  
改正了解析电模型中的Lplasma与LD表达式；考虑线圈阻抗的几何校正  
v200421 by陈鹏 fix了无法使用LZS的veff的bug，与LZS的解析电磁模型一致  
  
v200424 多自变量运行模式及后处理可用；改用2014Cazzador拟合表达式计算vst  
v200425 新的ELISE case benchmark，见result_of_main_v200418.md与eP  
典型参数下适用与结果分析 参数化分析  
复电导率虚部的影响  
  
v200514 不同种类碰撞频率绘图。与2014Cazzador结果一致。  
  
v200607  
以2019Raunera-fig7为benchmark，基本一致，见Result\benchmark_nu_with2019Raunera.png  
2018Jainb-fig.41中vstoc与气压负相关(计算式中与气压无关)，且值偏大。  
气压加入参数化。  
CHARLIE实现一一对应参数耦合计算  
  
v201010 v0.1 
201030 开始Git版本管理  