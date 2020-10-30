# history
赵鹏，左晨v190111
陈鹏v200420

# 主脚本
main.m 主脚本，修改输入参数和后处理代码，执行
## equivalent_medium_model_of_plasma部分 的子函数
k_eniz.m+IONIZATION.txt 计算电子中性粒子电离碰撞反应系数的函数
k_enp.m+ELASTIC.txt 计算电子中性粒子弹性碰撞反应系数的函数
## electric model部分 的子函数
check_Nagaoka.m Nagaoka.txt 查表得长冈系数，对理想螺线管电感做校正得实际电感

# 单独脚本
solve_stoc_eqns.m 单独脚本，用于求解stoc方程组
temp_test_code.m 单独脚本，测试子函数