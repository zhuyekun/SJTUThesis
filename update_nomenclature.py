
content = r"""% !TEX root = ../main.tex

\begin{nomenclature}
    \label{chap:symb}

    \begin{longtable}{rl}
      $d$ & 空间维度 \\
      $D$ & 空间区域 $D\subset \mathbb{R}^d$ \\
      $\partial D$ & $D$ 的边界 \\
      $\sS^{d-1}$ & $\mathbb{R}^d$ 中的单位球 \\
      $S_{d-1}$ & $\sS^{d-1}$ 的表面积 \\
      $\br$ & 位置向量 \\
      $\bOmega$ & 方向单位向量 \\
      $\bm{v}$ & 速度向量 $\bm{v}=\lvert \bm{v} \rvert \, \bOmega$ \\
      $I$ & 辐射强度 \\
      $I_{-}$ & $\Gamma_{-}$ 上的入流边界强度 \\
      $\mut$ & 总截面 \\
      $\mus$ & 散射截面 \\
      $\mu_a$ & 吸收截面 $\mu_a=\mut-\mus$ \\
      $p$ & 相函数 \\
      $\bn(\br)$ & $\br\in\partial D$ 处的单位外法向量 \\
      $\Gamma_{\pm}$ & 入流/出流边界集合 \\
      $s$ & 特征线上的路径长度参数 \\
      $s_{-}(\br,\bOmega)$ & 沿 $-\bOmega$ 方向从 $\br$ 到边界的距离 \\
      $g$ & Henyey--Greenstein 相函数的不对称参数 \\
      $\Phi$ & 标量密度（$I$ 的角度平均值） \\
      $\tau_{\br,\bOmega}(s_1, s_2)$ & 光学深度 \\
      $G$ & RTE 的格林函数 \\
      $\A$ & 解算子 $\A[I_-]=I$ \\
      $\J$ & 边界（衰减）算子 \\
      $\cL$ & 提升算子（带 $\mu_s$ 的特征积分） \\
      $\cS$ & 散射算子 \\
      $\rho_p$ & 加权 $L^p$ 空间中的谱半径 \\
      $\sigma$ & 特征长度 \\
      $\zeta_\sigma$ & 磨光核 \\
      $\delta_{\{\br'\}}$ & $\br'$ 处 $\partial D$ 上的狄拉克分布 \\
      $\delta$ & 方向空间上的狄拉克 $\delta$ 函数 \\
      $\br^{\text{mesh}}_i$ & 空间网格点 \\
      $(\mut^{\text{mesh}})_i$ & 空间依赖的总截面 \\
      $(\mus^{\text{mesh}})_i$ & 空间依赖的散射截面 \\
      $D_\mu$ & 截面子域 \\
      $d_{\text{model}}$ & 潜在/截断表示维度 \\
      $d_{\text{mlp}}$ & MLP 层的宽度（隐藏层大小） \\
      $N_{\text{mesh}}$ & 网格点数量 \\
      $N_{\text{quad}}$ & 角度求积点数量 \\
      $N_{\text{mlp}}$ & 衰减模块中的 MLP 层数 \\
      $N_{\ell}$ & 散射（残差）块的数量 \\
      $H$ & 注意力头数 \\
      $d_k, d_v$ & 注意力键/值嵌入维度 \\
      $d_{\tau}$ & 光学深度特征向量的维度 \\
      $\bG^{\text{NN}}$ & 格林函数的离散向量表示 \\
      $\bm{W}^{\ell}$ & 第 $\ell$ 个散射块中的权重矩阵 \\
      $\tau^{\text{NN}}_{-}$ & 神经光学深度估计 \\
      $\tau_{-,t}$ & 到边界的总光学深度 \\
      $\tau_{-,s}$ & 到边界的散射光学深度 \\
      $\bm{b}^{\ell}$ & 第 $\ell$ 个散射块中的偏置向量 \\
      $w_i$ & 求积权重 \\
      $\ell$ & 均方误差 \\
      $\mathcal{L}$ & 经验训练损失 \\
      $B$ & 批量大小 \\
      $\eta$ & 学习率 \\
    \end{longtable}

\end{nomenclature}
"""

with open('/root/projects/latex/SJTUThesis/contents/nomenclature.tex', 'w', encoding='utf-8') as f:
    f.write(content)
print("Successfully wrote contents/nomenclature.tex")
