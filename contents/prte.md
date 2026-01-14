# 辐射传输方程中伴随格林函数的定义、推导与应用
辐射传输方程（Radiative Transfer Equation, RTE）是描述辐射在介质中传输过程（包括衰减、散射与发射）的核心控制方程，广泛应用于大气科学、遥感探测、核工程等领域。对于含体源的辐射传输问题，格林函数法是求解非齐次方程的有效工具，其中伴随格林函数因能简化边界条件处理与反演问题求解，成为辐射传输理论中的关键概念。本文基于指定的符号系统，系统阐述辐射传输方程中伴随格林函数的定义、严格推导过程及典型应用场景。

## 1 基本方程与符号约定
考虑含体源项的稳态辐射传输方程，其标准形式为：
$$
(\boldsymbol{\Omega}\cdot\nabla + \mu_t) I(\boldsymbol{r},\boldsymbol{\Omega}) - \mu_s \int_{\mathbb{S}} p(\boldsymbol{\Omega},\boldsymbol{\Omega}') I(\boldsymbol{r},\boldsymbol{\Omega}')\,\mathrm{d}\boldsymbol{\Omega}' = q(\boldsymbol{r},\boldsymbol{\Omega}),
$$
其中各符号的物理意义与数学定义如下：
- $\boldsymbol{r}$ 为空间位置矢量；
- $\boldsymbol{\Omega}$ 为辐射传播方向单位矢量，取值于单位球面 $\mathbb{S}$（立体角域）；
- $I(\boldsymbol{r},\boldsymbol{\Omega})$ 为辐射强度，是位置与方向的二元函数；
- $\mu_t = \mu_a + \mu_s$ 为总衰减系数，$\mu_a$ 为吸收系数，$\mu_s$ 为散射系数；
- $p(\boldsymbol{\Omega},\boldsymbol{\Omega}')$ 为散射相函数，描述辐射从方向 $\boldsymbol{\Omega}'$ 散射至 $\boldsymbol{\Omega}$ 的概率，满足归一化条件 $\int_{\mathbb{S}} p(\boldsymbol{\Omega},\boldsymbol{\Omega}')\mathrm{d}\boldsymbol{\Omega}' = 1$；
- $q(\boldsymbol{r},\boldsymbol{\Omega})$ 为体源项，表征单位体积、单位立体角内的辐射发射强度。

为简化方程表述，定义线性辐射传输算子 $\mathcal{L}$，对任意辐射强度函数 $I(\boldsymbol{r},\boldsymbol{\Omega})$，算子作用为：
$$
\mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] = (\boldsymbol{\Omega}\cdot\nabla + \mu_t) I(\boldsymbol{r},\boldsymbol{\Omega}) - \mu_s \int_{\mathbb{S}} p(\boldsymbol{\Omega},\boldsymbol{\Omega}') I(\boldsymbol{r},\boldsymbol{\Omega}')\,\mathrm{d}\boldsymbol{\Omega}'.
$$
基于此，含体源的辐射传输方程可简洁表示为算子形式：
$$
\mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] = q(\boldsymbol{r},\boldsymbol{\Omega}).
$$
该方程为非齐次线性积分-微分方程，其求解的核心思路是利用格林函数的脉冲响应特性，将任意体源项 $q(\boldsymbol{r},\boldsymbol{\Omega})$ 分解为无数脉冲源的叠加，再通过线性算子的叠加性得到总响应。由于直接求解正算子 $\mathcal{L}$ 的格林函数会面临复杂的边界条件处理问题，实际中通常采用伴随算子与伴随格林函数的框架完成求解。

## 2 伴随格林函数的定义
伴随格林函数的定义建立在伴随算子与脉冲源响应的基础上，首先明确伴随算子的定义：对于任意两个满足边界条件的光滑函数 $I(\boldsymbol{r},\boldsymbol{\Omega})$ 与 $G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')$，伴随算子 $\mathcal{L}^*$ 满足内积恒等式：
$$
\int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \cdot \mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r} = \int_{\mathcal{V}} \int_{\mathbb{S}} I(\boldsymbol{r},\boldsymbol{\Omega}) \cdot \mathcal{L}^*[G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')] \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r} + \int_{\partial\mathcal{V}} \int_{\mathbb{S}} (\boldsymbol{\Omega}\cdot\boldsymbol{n}) I(\boldsymbol{r},\boldsymbol{\Omega}) G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}\mathrm{d}S,
$$
其中 $\mathcal{V}$ 为求解空间区域，$\partial\mathcal{V}$ 为区域边界，$\boldsymbol{n}$ 为边界外法向单位矢量，$\mathrm{d}S$ 为边界面积元，等式右侧的面积分即为边界项，源于微分算子的散度转化。

基于上述伴随算子定义，辐射传输方程的伴随格林函数 $G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')$ 被定义为满足以下脉冲源条件的函数：
$$
\mathcal{L}^*[G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')] = \delta(\boldsymbol{r} - \boldsymbol{r}') \delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}'),
$$
其中 $\delta(\boldsymbol{r} - \boldsymbol{r}')$ 为三维空间狄拉克δ函数，$\delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}')$ 为方向空间狄拉克δ函数，二者共同构成“位置-方向”二维脉冲源，表征仅在位置 $\boldsymbol{r}'$、方向 $\boldsymbol{\Omega}'$ 处存在单位强度的辐射源，其余位置与方向的源强度为零。

该式即为伴随格林函数的核心定义式，其物理意义可表述为：$G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')$ 是伴随算子 $\mathcal{L}^*$ 对位于 $(\boldsymbol{r}',\boldsymbol{\Omega}')$ 的单位脉冲源的响应，即当在 $(\boldsymbol{r}',\boldsymbol{\Omega}')$ 处施加单位脉冲辐射源时，在观测点 $(\boldsymbol{r},\boldsymbol{\Omega})$ 处产生的辐射强度。这一定义是辐射传输格林函数法的基础约定，广泛见于线性输运理论与辐射传输相关文献中。

结合前文正算子 $\mathcal{L}$ 的表达式，通过内积恒等式的展开可进一步得到伴随算子 $\mathcal{L}^*$ 的具体形式：
$$
\mathcal{L}^*[G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')] = -(\boldsymbol{\Omega}\cdot\nabla) G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') + \mu_t G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') - \mu_s \int_{\mathbb{S}} p(\boldsymbol{\Omega}',\boldsymbol{\Omega}) G(\boldsymbol{r},\boldsymbol{\Omega}';\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}'.
$$
对比正算子 $\mathcal{L}$ 与伴随算子 $\mathcal{L}^*$ 的表达式可见，二者的核心差异在于微分项的符号（正算子为 $+\boldsymbol{\Omega}\cdot\nabla$，伴随算子为 $-\boldsymbol{\Omega}\cdot\nabla$），这一差异源于方向梯度的散度转化过程，也使得伴随格林函数对应辐射的“反向传输”过程，从而大幅简化边界条件的处理。

## 3 伴随格林函数的推导与积分解的构建
伴随格林函数的推导核心是利用内积恒等式与δ函数的筛选性质，将非齐次方程 $\mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] = q(\boldsymbol{r},\boldsymbol{\Omega})$ 转化为积分形式，具体步骤如下：

第一步，对非齐次方程施加格林函数加权积分。将方程 $\mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] = q(\boldsymbol{r},\boldsymbol{\Omega})$ 两边乘以伴随格林函数 $G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')$，并在求解区域 $\mathcal{V}$ 与方向空间 $\mathbb{S}$ 上积分，得到：
$$
\int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \cdot \mathcal{L}[I(\boldsymbol{r},\boldsymbol{\Omega})] \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r} = \int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \cdot q(\boldsymbol{r},\boldsymbol{\Omega}) \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r}.
$$

第二步，代入伴随算子内积恒等式。将伴随算子的内积恒等式代入上式左侧，将左侧积分分解为伴随算子的内积项与边界项，即：
$$
\int_{\mathcal{V}} \int_{\mathbb{S}} I(\boldsymbol{r},\boldsymbol{\Omega}) \cdot \mathcal{L}^*[G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')] \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r} + \int_{\partial\mathcal{V}} \int_{\mathbb{S}} (\boldsymbol{\Omega}\cdot\boldsymbol{n}) I(\boldsymbol{r},\boldsymbol{\Omega}) G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}\mathrm{d}S = \int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \cdot q(\boldsymbol{r},\boldsymbol{\Omega}) \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r}.
$$

第三步，利用伴随格林函数的脉冲源条件化简积分。根据伴随格林函数的定义式 $\mathcal{L}^*[G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}')] = \delta(\boldsymbol{r} - \boldsymbol{r}') \delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}')$，将其代入上式左侧的体积分，利用狄拉克δ函数的筛选性质 $\int_{\mathcal{V}} f(\boldsymbol{r}) \delta(\boldsymbol{r} - \boldsymbol{r}') \mathrm{d}\boldsymbol{r} = f(\boldsymbol{r}')$ 与 $\int_{\mathbb{S}} f(\boldsymbol{\Omega}) \delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega} = f(\boldsymbol{\Omega}')$，可直接提取出待求的辐射强度 $I(\boldsymbol{r}',\boldsymbol{\Omega}')$，即：
$$
\int_{\mathcal{V}} \int_{\mathbb{S}} I(\boldsymbol{r},\boldsymbol{\Omega}) \cdot \delta(\boldsymbol{r} - \boldsymbol{r}') \delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}\mathrm{d}\boldsymbol{r} = I(\boldsymbol{r}',\boldsymbol{\Omega}').
$$

第四步，整理得到积分解。将上式代入积分恒等式，并将观测点符号 $(\boldsymbol{r}',\boldsymbol{\Omega}')$ 替换为常规表示 $(\boldsymbol{r},\boldsymbol{\Omega})$，源点符号保持 $(\boldsymbol{r}',\boldsymbol{\Omega}')$ 不变，最终得到辐射传输方程的格林函数积分解：
$$
I(\boldsymbol{r},\boldsymbol{\Omega}) = \int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') q(\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}'\mathrm{d}\boldsymbol{r}' - \int_{\partial\mathcal{V}} \int_{\mathbb{S}} (\boldsymbol{\Omega}'\cdot\boldsymbol{n}') I(\boldsymbol{r}',\boldsymbol{\Omega}') G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}'\mathrm{d}S'.
$$

上述推导过程中，散度定理的应用是实现微分算子向边界积分转化的关键。对于正算子中方向微分项 $\boldsymbol{\Omega}\cdot\nabla I$ 的处理，通过矢量恒等式 $\boldsymbol{\Omega}\cdot\nabla(GI) = G(\boldsymbol{\Omega}\cdot\nabla I) + I(\boldsymbol{\Omega}\cdot\nabla G)$ 进行变形，再对 $\boldsymbol{\Omega}\cdot\nabla(GI)$ 应用散度定理 $\int_{\mathcal{V}} \boldsymbol{\Omega}\cdot\nabla(GI) \mathrm{d}\boldsymbol{r} = \int_{\partial\mathcal{V}} (\boldsymbol{\Omega}\cdot\boldsymbol{n}) GI \mathrm{d}S$，完成体积分向边界积分的转化，这一过程也是伴随算子边界项的唯一来源。

该积分解中包含两项核心贡献：第一项为体源项 $q(\boldsymbol{r}',\boldsymbol{\Omega}')$ 对应的体积分贡献，体现了任意体源可分解为无数脉冲源的叠加，总响应为各脉冲源响应（伴随格林函数）的积分叠加；第二项为边界项，表征边界上的辐射强度对区域内辐射场的贡献，其符号由方向与边界法向的夹角决定，$\boldsymbol{\Omega}'\cdot\boldsymbol{n}' > 0$ 对应辐射离开边界，$\boldsymbol{\Omega}'\cdot\boldsymbol{n}' < 0$ 对应辐射入射至边界。

在无边界源（如真空边界或黑体边界无入射辐射）的情况下，边界项为零，积分解简化为仅含体积分的形式：
$$
I(\boldsymbol{r},\boldsymbol{\Omega}) = \int_{\mathcal{V}} \int_{\mathbb{S}} G(\boldsymbol{r},\boldsymbol{\Omega};\boldsymbol{r}',\boldsymbol{\Omega}') q(\boldsymbol{r}',\boldsymbol{\Omega}') \mathrm{d}\boldsymbol{\Omega}'\mathrm{d}\boldsymbol{r}'.
$$

## 4 伴随格林函数的应用
伴随格林函数凭借其对边界条件的简化处理能力与脉冲响应的叠加特性，在辐射传输的正演计算与反演问题中具有广泛应用，主要体现在以下三个核心领域：

### 4.1 复杂介质辐射传输的正演计算
在三维非均匀介质（如大气气溶胶、植被冠层、非均匀水体）的辐射传输问题中，解析求解辐射传输方程通常较为困难，伴随格林函数法可将复杂的积分-微分方程转化为积分形式，通过数值积分或蒙特卡洛模拟求解。例如，在大气遥感领域，利用伴随格林函数可高效计算非均匀大气中体源（如大气分子发射、气溶胶散射源）的辐射场分布，为卫星观测数据的模拟提供理论基础。对于时变或频域辐射传输问题，伴随格林函数还可通过傅里叶变换扩展至相应域，实现动态辐射过程的正演模拟。

### 4.2 辐射传输反演问题的求解
反演问题（如从辐射观测数据反演介质的光学参数 $\mu_t$、$\mu_s$ 或体源分布 $q$）是辐射传输理论的重要应用方向，伴随格林函数是构建反演算法的关键工具。由于反演问题通常为病态问题，需要利用梯度信息进行优化求解，而伴随格林函数可直接给出观测值对反演参数的灵敏度（梯度）。例如，在地表遥感反演中，通过伴随格林函数可高效计算地表反射率、植被覆盖率等参数对卫星观测辐射强度的影响，显著降低反演算法的计算复杂度。

### 4.3 边界条件敏感问题的分析
在核辐射防护、高温燃烧等领域，边界条件（如壁面反射、入射辐射强度）对辐射场分布具有显著影响。伴随格林函数对应的反向传输特性，可直观反映边界辐射对区域内任意观测点的贡献，从而实现边界条件敏感性的定量分析。例如，在核反应堆的辐射屏蔽设计中，利用伴随格林函数可计算反应堆壁面泄漏辐射对周围环境的辐射剂量分布，为屏蔽结构的优化提供理论依据。

此外，伴随格林函数还可与球谐函数展开、有限体积法等数值方法结合，形成高效的辐射传输数值求解框架。例如，通过球谐函数展开伴随格林函数，可将积分方程转化为线性代数方程组，大幅提升求解效率；在三维复杂几何区域中，结合有限体积法离散伴随格林函数的积分形式，可实现复杂边界条件下的高精度辐射场计算。

## 5 结论
辐射传输方程的伴随格林函数是基于伴随算子定义的脉冲响应函数，其核心定义式 $\mathcal{L}^*[G] = \delta(\boldsymbol{r} - \boldsymbol{r}') \delta(\boldsymbol{\Omega} - \boldsymbol{\Omega}')$ 源于线性算子的脉冲响应叠加思想，通过将任意体源分解为无数脉冲源的叠加，利用线性算子的叠加性得到总响应。伴随格林函数的推导过程借助伴随算子的内积恒等式与散度定理，将非齐次辐射传输方程转化为积分形式，其中δ函数的筛选性质是提取待求辐射强度的关键，而散度定理则完成了微分算子向边界积分的转化。

在应用层面，伴随格林函数不仅简化了复杂介质辐射传输的正演计算，更为辐射传输反演问题与边界条件敏感分析提供了高效的理论工具，是辐射传输理论与工程应用之间的重要桥梁。后续研究中，伴随格林函数与机器学习、高性能计算的结合，将进一步拓展其在复杂多介质、多物理场耦合辐射传输问题中的应用范围。

## 参考文献
[1] Case, K. M., & Zweifel, P. F. (1967). *Linear Transport Theory*. Addison-Wesley.

[2] Lyapustin, A., & Knyazikhin, Y. (2002). Green's function method in the radiative transfer problem. II. Spatially heterogeneous anisotropic surface. *Applied Optics*, 41(27), 5600-5606. https://doi.org/10.1364/AO.41.005600

[3] Duderstadt, J. J., & Martin, W. R. (1979). *Transport Theory*. John Wiley & Sons.

[4] Modest, M. F. (2021). *Radiative Heat Transfer* (4th ed.). Academic Press/Elsevier.

---
需要我帮你将这份文档转换成**期刊投稿格式的tex文件**，或者补充**数值算例的公式推导**吗？

关键点：对固定的方向 $\bOmega\in\mathbb S^{d-1}$，$\bOmega$ 与空间变量 $\br$ 无关，是“常向量”。令 $f(\br)=G(\br,\bOmega;\br',\bOmega')\,I(\br,\bOmega)$，则有恒等式
[
\nabla\cdot(\bOmega f)=\bOmega\cdot\nabla f + f,(\nabla\cdot\bOmega).
]
但因为 $\bOmega$ 不依赖于 $\br$，所以 $\nabla\cdot\bOmega=0$，于是
[
\nabla\cdot(\bOmega f)=\bOmega\cdot\nabla f.
]
写成分量更直观：设 $\br=(x_1,\dots,x_d)$，$\bOmega=(\Omega_1,\dots,\Omega_d)$，
[
\nabla\cdot(\bOmega f)=\sum_{i=1}^d \frac{\partial}{\partial x_i}(\Omega_i f)
=\sum_{i=1}^d \Omega_i \frac{\partial f}{\partial x_i}
=\bOmega\cdot\nabla f,
]