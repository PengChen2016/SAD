# SAD  
**the Simple Analysis of the Driver in RF hydrogen ion source**  

The code is still under development, though the the existing modules have been compared with the results reported by the literature. (See the unit test cases, and the results are listed in .\others\Figures during code developing.pptx in Chinese. Sorry for nonstandard references.) 

If this repository is helpful for you, please email PengChen(pengchen.hust@outlook.com) just let me know. It is my honor to know that.  
Please feel free to use the Discussions and pull request in Github.  
The data and the related document(in Chinese) are available from the author upon reasonable request. 

## Content
 ### Examples of the plasma models and other electric models
.\CHARLIE_sweep210415.m
.\BUG_sweep210508.m
.\CHARLIE_nonuniform210426.m
 
#### results and analysis
.\others\for_paper210415\
 
### Examples of the FEM model
.\others\for_paper210415\BUG1.aedt
.\others\for_paper210415\CHARLIE_1MHz.aedt 

## Reference
[1]	Zielke D P, Briefi S, and Fantz U 2021 J. Phys. D: Appl. Phys. 54 155202, https://doi.org/10.1088/1361-6463/abd8ee
[2]	Kralkina E A, Rukhadze A A, Pavlov V B, Vavilin K V, Nekliudova P A, Petrov A K and Alexandrov A F 2016 Plasma Sources Sci. Technol. 25 015016
[3]	Kraus W, Fantz U, Heinemann B and Franzen P. 2015 Fusion Eng. Des. 91 16, https://doi.org/10.1016/j.fusengdes.2014.11.015
[4]	Rauner D, Briefi S and Fantz U 2017 Plasma Sources Sci. Technol. 26 095004
[5]	Rauner D, Briefi S and Fantz U. 2019 Plasma Sources Sci. Technol. 28 095011
[6]	Piejak R B, Godyak V A and Alexandrovich B M 1992 Plasma Sources Sci. Technol. 1 179
[7]	Cazzador M 2014 MSc Thesis University of Padua, Padua
[8]	Nishida K, Mochizuki S, Ohta M, Yasumoto M, Lettry J, Mattei S, Mattei S and Hatayama A 2014 Rev. Sci. Instrum. 85 02B117
[9]	Sudhir D, Bandyopadhyay M, Kraus W, Gahlaut A, Bansal G and Chakraborty A 2014 Rev. Sci. Instrum. 85 013510
[10]	Bandyopadhyay M, Sudhir D and Chakraborty A 2015 Nucl. Fusion 55 033017
[11]	Jain P, Recchia M, Veltri P, Cavenago M, Maistrello A and Gaio E 2018 IEEE Access 6 29665
[12]	Chabert P and Braithwaite N 2011 Physics of Radio-Frequency Plasmas (Cambridge: Cambridge University Press)
[13]	Zhang L G, Chen D Z, Li D, Liu K F, Li X F, Pan R M and Fan M W Fusion Eng. Des. 2016 103 74
[14]	Lee S W, Goulding R H, Kang Y W, Shin K and Welton R F 2010 Rev. Sci. Instrum. 81 02A726
[15]	Grudiev A, Lettry J, Mattei S, Paoluzzi M and Scrivens R Rev. Sci. Instrum. 2014 85 02B134
[16]	Zhao P, Li D, Chen D Z, Zuo C, Li X F and Fan M W 2018 Fusion Eng. Des. 132 29
[17]	Chen P, Li D, Chen D Z, Song F, Zuo C, Zhao P, Lei G J 2018 AIP Conf Proc. 2052 040018
[18]	Zielke D P, Rauner D, Briefi S, Lishev S and Fantz U 2021 Plasma Sources Sci. Technol. in press, https://doi.org/10.1088/1361-6595/ac0396
[19]	Kobayashi W, Nishida K, Mattei S, Lettry J, Hoshino K and Hatayama A 2018 AIP Conf Proc. 2052 050009
[20]	Rauner D, Mattei S, Briefi S, Fantz U, Hatayama A, Lettry J, Nishida K and Tran M Q 2017 AIP Conf Proc. 1869 030035
[21]	Lieberman M A and Lichtenberg A J 2005 Principles of Plasma Discharges and Materials Processing 2nd edn (Hoboken, New Jersey: Wiley)
[22]	Yoon N S, Kim S S, Chang C S, Choi D I 1996 Phys. Rev. E 54 757
[23]	Vahedi V et al 1995 J. Appl. Phys. 78 1446
[24]	Nagaoka H. 1909 J. Coll. Sci., Imp. Univ. Tokyo 27 1, https://doi.org/10.15083/00037799
[25]	Phelps database, www.lxcat.net, retrieved on December 14, 2020.
[26]	Buckman S J and Phelps A V, 1995 JILA Information Center Report No. 27, University of Colorado
[27]	Rauner D 2019 PhD Thesis Augsburg University, Augsburg
[28]	Briefi S and Fantz U. 2018 AIP Conf Proc. 2052 040005, https://doi.org/10.1063/1.5083739
[29]	Mcneely P, Wunderlich W and the NNBI Team 2011 Plasma Sources Sci. Technol. 20 045005, http://doi.org/10.1088/0963-0252/20/4/045005
[30]	McNeely P, Dudin S V, Christ-Koch S and Fantz U. 2009 Plasma Sources Sci. Technol. 18 014011
[31]	Zhao P, Li D, Chen D, Liu K, Li X, Wang H, Zuo C and Fan M. 2017 Nucl. Instrum. Methods Phys. Res., Sect. A 858 90, http://dx.doi.org/10.1016/j.nima.2017.03.054


## get collision cross section
https://nl.lxcat.net/data/set_type.php -> type: Scattering crossing sections -> database: all -> specA: Ground states: e -> specB: Ground states: Ar -> groups: Elastic, Ionization -> processes: Biagi-v7.1, Morgan -> compare and download -> put each data set into a single file, record the reference at the first line -> modify the get_Xsec.m  
Cross sections of e-H2, e-Ar has been included in the ./packages/physic_tools.  

## compare data
output data by get_output_json(), then compare file by arbitrary text comparator, such as plug-in of code editor

## History before version control by Git (in Chinese)
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