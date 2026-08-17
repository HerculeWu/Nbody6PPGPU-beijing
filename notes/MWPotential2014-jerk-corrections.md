# MWPotential2014 的两处 jerk 公式修正

## 1. 修正范围与物理背景

Nbody6++ Beijing 在 `KZ(14)=5` 时使用与 Galpy
`MWPotential2014` 相同的三分量银河系势：

1. 带指数截断的球对称幂律核球；
2. Miyamoto–Nagai 盘；
3. NFW 暗物质晕。

本次修正不改变势函数 \(\Phi\)、加速度
\(\boldsymbol{a}=-\nabla\Phi\)、各分量参数或在
\(R_0=8\,\mathrm{kpc}\)、\(V_0=220\,\mathrm{km\,s^{-1}}\)
处的归一化。修正对象仅为加速度的一阶时间导数

\[
\boldsymbol{j}\equiv\frac{\mathrm d\boldsymbol{a}}{\mathrm dt},
\]

即 Hermite 积分器使用的 jerk。

受影响的调用路径为：

```text
F_PowCut / F_NFW
        ↓
      fmwpot
        ↓
XTRNLF（恒星相对质心的差分外场）
GCINT （星团质心的银河轨道）
XTRNLT（全局外场调用）
```

因此，旧代码描述的静态 MWPotential2014 势和加速度是正确的，
但其时间导数与该加速度场不完全相容。

## 2. 球对称场的通用 jerk 公式

对任意球对称势，可将加速度写为

\[
\boldsymbol{a}(\boldsymbol{r})=A(r)\boldsymbol{r},
\qquad
r=\lvert\boldsymbol{r}\rvert .
\]

令

\[
u=\boldsymbol{r}\cdot\boldsymbol{v},
\qquad
\dot r=\frac{u}{r}.
\]

沿粒子轨道对加速度求导：

\[
\begin{aligned}
\boldsymbol{j}
 &=\frac{\mathrm d}{\mathrm dt}
   \left[A(r)\boldsymbol{r}\right] \\
 &=A(r)\boldsymbol{v}
   +A'(r)\dot r\,\boldsymbol{r} \\
 &=A(r)\boldsymbol{v}
   +\frac{A'(r)}{r}
    (\boldsymbol{r}\cdot\boldsymbol{v})\boldsymbol{r}.
\end{aligned}
\]

这条关系是下面两处修正的共同出发点。它也说明 jerk
不仅包含 \(A\boldsymbol{v}\)，还必须包含粒子径向运动穿过
空间变化外场时产生的 \(A'(r)\dot r\boldsymbol{r}\)。

## 3. 错误一：截断幂律核球的密度梯度项多除以 \(r^2\)

### 3.1 文件位置

- 文件：`src/Main/fpowercut.f`
- 修正后行号：第 27 行
- 子程序：`F_PowCut`

### 3.2 势模型

核球密度为

\[
\rho(r)=\rho_0 r^{-\alpha}
        \exp\left[-\left(\frac{r}{r_c}\right)^2\right],
\]

包围质量为

\[
M(r)=4\pi\int_0^r \rho(s)s^2\,\mathrm ds.
\]

球对称引力加速度写成

\[
\boldsymbol{a}_{\mathrm b}
 =-\frac{M(r)}{r^3}\boldsymbol{r}.
\]

对应代码变量为

\[
\texttt{mr}=-\frac{M(r)}{r^3},
\qquad
\texttt{rvr2}
 =q=\frac{\boldsymbol{r}\cdot\boldsymbol{v}}{r^2},
\qquad
\texttt{dmr}=-4\pi\rho(r).
\]

### 3.3 数学推导

由

\[
\frac{\mathrm dM}{\mathrm dr}=4\pi r^2\rho(r),
\qquad
\dot r=qr,
\]

可得

\[
\dot M=4\pi r^3\rho(r)q.
\]

对 \(A(r)=-M(r)r^{-3}\) 求时间导数：

\[
\begin{aligned}
\dot A
 &=-\dot M r^{-3}
   +3M r^{-4}\dot r \\
 &=-4\pi\rho(r)q
   +3\frac{M(r)}{r^3}q \\
 &=-4\pi\rho(r)q-3Aq.
\end{aligned}
\]

代回通用 jerk 公式：

\[
\boxed{
\boldsymbol{j}_{\mathrm b}
=A\left(\boldsymbol{v}-3q\boldsymbol{r}\right)
-4\pi\rho(r)q\boldsymbol{r}
}.
\]

按代码变量写为

```fortran
FD = mr*(XIDOT - 3*rvr2*XI) + dmr*rvr2*XI
```

### 3.4 原错误及其物理含义

旧代码为

```fortran
FD = mr*(XIDOT - 3*rvr2*XI) + dmr*rvr2/r2*XI
```

其中密度项被额外除以 \(r^2\)。但
\(\mathrm dM/\mathrm dr=4\pi r^2\rho\) 中的 \(r^2\)
已经与球对称力中的 \(r^{-3}\) 及
\(\dot r=qr\) 完全约掉，不应再次出现 \(r^{-2}\)。

这个额外因子会人为改变核球 jerk 的径向尺度：

- 在 \(r<1\) 的代码长度单位内放大密度梯度贡献；
- 在 \(r>1\) 时压低该贡献；
- 使 jerk 不再等于同一加速度场沿真实轨道的时间导数。

### 3.5 实际修复

删除密度项中的 `/r2`：

```diff
- FD = mr*(XIDOT -3.0*rvr2*XI) + dmr*rvr2/r2*XI
+ FD = mr*(XIDOT -3.0*rvr2*XI) + dmr*rvr2*XI
```

## 4. 错误二：NFW 径向导数的系数和组合符号错误

### 4.1 文件位置

- 文件：`src/Main/fnfw.f`
- 修正后行号：第 23–24 行和第 28 行
- 子程序：`F_NFW`

### 4.2 势和加速度

NFW 势写为

\[
\Phi_{\mathrm h}(r)
=-H\frac{\ln(1+r/a)}{r},
\]

其中 \(H=\texttt{H\_AMP}\)，\(a=\texttt{H\_A}\)。
定义

\[
L=\ln(1+r/a),
\qquad
s=a+r.
\]

其加速度可写成

\[
\boldsymbol{a}_{\mathrm h}=A(r)\boldsymbol{r},
\]

其中

\[
A(r)
=H\left[
\frac{1}{r^2s}
-\frac{L}{r^3}
\right].
\]

这正是代码变量 `DP`：

```fortran
DP = H_AMP*(1/(r2*ra) - log(1+r/H_A)/r3)
```

### 4.3 数学推导

分别求导：

\[
\frac{\mathrm d}{\mathrm dr}
\left(\frac{1}{r^2s}\right)
=-\frac{2}{r^3s}-\frac{1}{r^2s^2},
\]

以及

\[
\frac{\mathrm d}{\mathrm dr}
\left(-\frac{L}{r^3}\right)
=-\frac{1}{r^3s}+\frac{3L}{r^4}.
\]

两项合并后，

\[
A'(r)
=H\left[
-\frac{3}{r^3s}
-\frac{1}{r^2s^2}
+\frac{3L}{r^4}
\right].
\]

令 \(u=\boldsymbol{r}\cdot\boldsymbol{v}\)，则

\[
\begin{aligned}
\dot A
 &=A'(r)\frac{u}{r} \\
 &=H\left[
-\frac{3}{r^4s}
-\frac{1}{r^3s^2}
+\frac{3L}{r^5}
\right]u.
\end{aligned}
\]

按照原代码便于实现的因式分解形式：

\[
\boxed{
\dot A
=H\left[
-\frac{3}{rs}
-\frac{1}{s^2}
+\frac{3L}{r^2}
\right]\frac{u}{r^3}
}.
\]

因此正确 jerk 为

\[
\boxed{
\boldsymbol{j}_{\mathrm h}
=A\boldsymbol{v}+\dot A\boldsymbol{r}
}.
\]

### 4.4 原错误及其物理含义

旧代码有两个独立问题：

1. `DDP` 第一项使用 `-2/(r*ra)`，遗漏了
   \(L'(r)=1/(a+r)\) 对 \(-L/r^3\) 的贡献；
   合并后的正确系数应为 `-3`。
2. `DDP` 按其定义就是 \(\dot A\)，所以 jerk 中应当
   加上 `DDP*XI`；旧代码却将其减去。

旧实现因此不满足

\[
\boldsymbol{j}
=(\boldsymbol{v}\cdot\nabla)\boldsymbol{a}
\]

这一静态外场必须满足的运动学关系。误差只作用于
与径向运动相关的部分：当
\(\boldsymbol{r}\cdot\boldsymbol{v}=0\) 时，
\(\dot A=0\)，所以严格圆轨道可能掩盖该问题；偏心轨道、
径向迁移以及星团内部恒星相对银河中心的径向速度都会使其显现。

### 4.5 实际修复

将径向导数系数从 `-2` 改为 `-3`，并按乘积求导结果
把 `DDP*XI` 加入 jerk：

```diff
- DDP = H_AMP*(-2/(r*ra) - 1/(ra*ra) + 3*log(1+r/H_A)/r2)
+ DDP = H_AMP*(-3/(r*ra) - 1/(ra*ra) + 3*log(1+r/H_A)/r2)

- FD = DP*XIDOT - DDP*XI
+ FD = DP*XIDOT + DDP*XI
```

## 5. 单位缩放与本次修正的关系

`fmwpot` 将 N-body 坐标转换为：

- 位置：kpc；
- 速度：pc/Myr；
- 加速度：pc/Myr²。

由于位置和速度在送入分量函数时使用了相差 \(10^3\)
的长度单位，分量函数返回的原始 jerk 也带有相应的
\(10^3\) 数值因子；`FD_PcMyr3` 中已有同样的
\(10^3\) 因子负责在返回 N-body 单位时抵消它。

本次错误不是这个既有单位换算造成的，也不需要改变该换算。
修正只恢复 \(\mathrm d\boldsymbol{a}/\mathrm dt\)
对空间和速度的正确解析依赖。

## 6. 对动力学的物理影响

Nbody6++ 的四阶 Hermite 预测—校正依赖同时一致的
\(\boldsymbol{a}\) 和 \(\boldsymbol{j}\)。旧公式不会改变
某一给定位置处的静态银河加速度，但会使预测出的下一时刻
外场偏离该静态势真实产生的变化率，可能表现为：

- 星团偏心银河轨道上的质心预测和校正误差增大；
- 恒星所受差分潮汐 jerk
  \(\boldsymbol{j}(\boldsymbol{R}+\boldsymbol{x})
   -\boldsymbol{j}(\boldsymbol{R})\) 不准确；
- 外场相关的轨道相位和长期能量行为与使用同一势函数、
  但直接重算加速度的积分方法出现额外差异。

修正后，核球和 NFW 分量的 jerk 都重新成为各自
MWPotential2014 加速度的解析时间导数。Miyamoto–Nagai
盘分量的原 jerk 推导是自洽的，因此未作修改。

