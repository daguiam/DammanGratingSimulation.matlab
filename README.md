# Damman_grating_simulation.matlab
点阵 [达曼光栅](https://www.researchgate.net/publication/294698009_Dammann_gratings_for_laser_beam_shaping) Matlab辅助设计

Author: [Ke Ding](https://github.com/Lonelyearner), *School of Optics and Photonics, BIT*

## 问题描述
寻找一种达曼光栅周期结构，使其在规定误差函数意义下最优
## 模型假设
光学系统是2f 系统，仿真时可以直接使用fft记算衍射图样
## 模型简化
- 周期结构参数(突变点对)取值在(0,1)
- 对周期结构密集采样，该过程造成的频谱损失忽略不计
- 达曼光栅相位为{0,𝜋}二值分布，对应于复振幅透射系数{1,−1}二值分布
- 二维达曼光栅可由一维达曼光栅生成，仅考虑一维光栅的优化
## 模型描述
- 变量初始化
- 单次迭代过程
   - 随机选取突变点对或扰动突变点
   - 对该光栅结构密集采样，采样点个数决定结构精度
   - 计算感兴趣衍射级次的衍射效率
   - 计算代价
- 使用模拟退火算法优化，获得一维达曼光栅结构最优解
- 生成二维达曼光栅，并显示结果
## 结果
