# H-DAB Dominance Extractor 工程文档

## 1. 项目概述

本项目实现了一个 Fiji/ImageJ 插件：

```text
H-DAB Dominance Extractor
```

插件安装位置：

```text
/Users/tom/Fiji/plugins/H_DAB_Dominance_Extractor.jar
```

Fiji 菜单入口：

```text
Plugins > Analyze > H-DAB Dominance Extractor
```

插件面向 H-DAB 免疫组化图像，目标不是替代 Fiji 自带的 Colour Deconvolution，而是在其基础思想上增加：

```text
1. 三张彩色 color deconvolution 图
2. 三张定量浓度图
3. DAB-positive mask
4. DAB dominance 筛选
5. 自动参数优化
6. 测量 CSV 输出
```

## 2. 需求定义

用户的原始需求可以拆成四个工程目标：

```text
1. 将 H-DAB 图像拆成 Hematoxylin、DAB、Residual 三张图
2. 输出的三张图应是 color deconvolution 风格的彩色图像
3. 自动区分哪些区域是真正 DAB 阳性
4. 对三张真实图的结果要优于传统算法
```

当前实现状态：

```text
已实现：彩色 H / DAB / Residual 输出
已实现：浓度图输出
已实现：DAB-positive mask
已实现：测量 CSV
已实现：三张真实图无标注 QC 对比
已实现：ImageJ 插件安装
```

## 3. 与 Fiji 自带 Colour Deconvolution 的关系

Fiji 自带插件：

```text
Image > Color > Colour Deconvolution
```

两者关系如下：

| 项目 | Fiji 自带 Colour Deconvolution | 本插件 |
|---|---|---|
| 核心分色模型 | Ruifrok-Johnston OD-space deconvolution | 同类 OD-space H-DAB deconvolution |
| H/DAB/Residual 输出 | 有 | 有 |
| 彩色 deconvolution 图 | 主要由自带插件显示方式决定 | 明确输出 `*_color.tif` |
| 浓度图 | 有类似 channel | 明确输出 float concentration TIFF |
| DAB 阳性判断 | 无 | 有 |
| DAB-Hematoxylin dominance | 无 | 有 |
| 自动阈值优化 | 无 | 有 |
| 自动测量 CSV | 无 | 有 |
| 目标 | 通用分色 | H-DAB DAB 阳性识别和定量 |

结论：

```text
分色底层原理相近；
本插件的工程价值主要在 DAB 阳性 mask、自动优化和定量输出。
```

## 4. 文件结构

核心源码：

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_src/H_DAB_Dominance_Extractor.java
```

测试入口：

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_src/H_DAB_Plugin_SmokeTest.java
```

插件菜单配置：

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_src/plugins.config
```

已安装插件：

```text
/Users/tom/Fiji/plugins/H_DAB_Dominance_Extractor.jar
```

真实图输出：

```text
/Users/tom/Fiji/cloud_extract_eval/real_plugin_outputs
```

真实图比较报告：

```text
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_report.md
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_qc_summary.csv
/Users/tom/Fiji/cloud_extract_eval/real_image_results/real_image_qc_metrics.csv
```

## 5. 输入输出

### 5.1 输入

输入图像要求：

```text
RGB image
Hematoxylin + DAB staining
JPG/TIF/PNG 均可，只要 ImageJ 能打开
```

已验证输入：

```text
/Users/tom/Desktop/NT_HTXR_K2_J1.jpg
/Users/tom/Desktop/slide-HFLday8-IRmodel-2026-03-26T12-37-47-R1-S10_11.3x.jpg
/Users/tom/Desktop/CB1_PCOS_Carotis_132_1_reexport.jpg
```

### 5.2 输出

每张图会输出：

```text
*_hematoxylin_color.tif
*_dab_color.tif
*_residual_color.tif

*_hematoxylin.tif
*_dab.tif
*_residual.tif

*_dab_dominance.tif
*_dab_positive_mask.tif
*_measurements.csv
```

输出解释：

| 文件 | 类型 | 用途 |
|---|---|---|
| `*_hematoxylin_color.tif` | RGB | 彩色 Hematoxylin-only deconvolution 图 |
| `*_dab_color.tif` | RGB | 彩色 DAB-only deconvolution 图 |
| `*_residual_color.tif` | RGB | 彩色 residual-only deconvolution 图 |
| `*_hematoxylin.tif` | float | Hematoxylin 浓度图，适合定量 |
| `*_dab.tif` | float | DAB 浓度图，适合定量 |
| `*_residual.tif` | float | residual 浓度图，适合 QC |
| `*_dab_dominance.tif` | float | `DAB - alpha * H` 优势分数 |
| `*_dab_positive_mask.tif` | 8-bit mask | DAB 阳性区域 |
| `*_measurements.csv` | CSV | 面积、IOD、平均 DAB 等定量结果 |

## 6. 算法流程

### 6.1 RGB 转 Optical Density

对每个像素：

```text
OD_R = -ln((R + 1) / 256)
OD_G = -ln((G + 1) / 256)
OD_B = -ln((B + 1) / 256)
```

`+1` 用于避免 `ln(0)`。

### 6.2 H-DAB stain vectors

默认 stain vectors：

```text
Hematoxylin = (0.650, 0.704, 0.286)
DAB         = (0.268, 0.570, 0.776)
```

这组向量来自常用 H-DAB color deconvolution preset。

Residual vector 由 H 和 DAB 向量推导：

```text
Residual_i = sqrt(max(0, 1 - H_i^2 - DAB_i^2))
```

然后归一化。

### 6.3 Stain matrix inversion

构造矩阵：

```text
M = [Hematoxylin, DAB, Residual]
```

对每个像素：

```text
[H, DAB, Residual] = inverse(M) * OD
```

负值裁剪为 0。

### 6.4 彩色 color deconvolution 图生成

每个 stain-only 彩色图通过该 stain 的浓度和 stain vector 重新投影回 RGB：

```text
R = 255 * exp(-amount * vector_R)
G = 255 * exp(-amount * vector_G)
B = 255 * exp(-amount * vector_B)
```

因此：

```text
Hematoxylin color = 只用 Hematoxylin 浓度重建
DAB color         = 只用 DAB 浓度重建
Residual color    = 只用 Residual 浓度重建
```

### 6.5 DAB dominance 阳性筛选

先按 tissue 区域内第 99 百分位归一化：

```text
H_norm
DAB_norm
Residual_norm
```

计算 DAB 优势分数：

```text
DAB_dominance = DAB_norm - alpha * H_norm
```

阳性 mask 条件：

```text
tissue
AND DAB_norm >= optimized_DAB_threshold
AND DAB_dominance >= optimized_dominance_threshold
AND DAB_norm >= 1.15 * H_norm
AND Residual_norm <= residual_threshold
```

### 6.6 高可信 DAB cloud floor

三张真实图里有大量淡褐色区域。为了避免把全部淡褐色背景都算成阳性，插件增加了一层高可信 DAB floor：

```text
DAB-dominant pixels = tissue AND DAB_norm >= 1.15 * H_norm

optimized_DAB_threshold 至少为：
85th percentile of DAB-dominant pixels
```

这个策略让插件更偏高可信 DAB 阳性区域，而不是广泛弱染色区域。

## 7. 参数

主要参数：

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `alpha` | 自动优化 | `DAB - alpha * H` 中的 H 惩罚权重 |
| `tissueMinOD` | 0.08 | tissue mask 的最低 OD |
| `dabMin` | 自动优化 | DAB 最低阈值 |
| `dominanceMin` | 自动优化 | DAB dominance 最低阈值 |
| `ratioMin` | 1.15 | DAB/Hematoxylin 最低比例 |
| `residualMax` | 1.05 | residual 上限 |
| `autoOptimize` | true | 是否自动优化参数 |

## 8. 安装与构建

### 8.1 构建命令

使用 Fiji 自带 JDK：

```text
/Users/tom/Fiji/java/macos-arm64/zulu21.42.19-ca-jdk21.0.7-macosx_aarch64/zulu-21.jdk/Contents/Home
```

编译和打包：

```bash
cd /Users/tom/Fiji
JDK="$PWD/java/macos-arm64/zulu21.42.19-ca-jdk21.0.7-macosx_aarch64/zulu-21.jdk/Contents/Home"
"$JDK/bin/javac" --release 8 -cp jars/ij-1.54p.jar \
  -d cloud_extract_eval/plugin_build/classes \
  cloud_extract_eval/plugin_src/H_DAB_Dominance_Extractor.java

cp cloud_extract_eval/plugin_src/plugins.config cloud_extract_eval/plugin_build/classes/plugins.config

"$JDK/bin/jar" cf plugins/H_DAB_Dominance_Extractor.jar \
  -C cloud_extract_eval/plugin_build/classes .
```

### 8.2 安装验证

JAR 内容应包含：

```text
H_DAB_Dominance_Extractor.class
plugins.config
```

菜单配置：

```text
Plugins>Analyze, "H-DAB Dominance Extractor", H_DAB_Dominance_Extractor
```

## 9. 测试与验证

### 9.1 Smoke test

测试入口：

```text
/Users/tom/Fiji/cloud_extract_eval/plugin_src/H_DAB_Plugin_SmokeTest.java
```

测试输出：

```text
/Users/tom/Fiji/cloud_extract_eval/color_deconv_test
```

已确认生成：

```text
CB1_color_deconv_hematoxylin_color.tif
CB1_color_deconv_dab_color.tif
CB1_color_deconv_residual_color.tif
CB1_color_deconv_hematoxylin.tif
CB1_color_deconv_dab.tif
CB1_color_deconv_residual.tif
CB1_color_deconv_dab_dominance.tif
CB1_color_deconv_dab_positive_mask.tif
CB1_color_deconv_measurements.csv
```

### 9.2 三张真实图验证

真实图插件输出：

```text
/Users/tom/Fiji/cloud_extract_eval/real_plugin_outputs
```

每张真实图已生成 9 个输出文件：

```text
3 color deconvolution TIFF
3 concentration TIFF
dominance TIFF
DAB-positive mask TIFF
measurement CSV
```

### 9.3 与传统方法比较

比较脚本：

```text
/Users/tom/Fiji/cloud_extract_eval/real_image_comparison.py
```

比较对象：

```text
Proposed H-DAB dominance
CCE-style auto seeds
Fixed H-DAB deconvolution
OD k-means
RGB k-means
```

真实三图无标注 QC 结果：

| Method | QC score | DAB area fraction | White positive % | Artifact positive % | DAB pos/non-pos ratio |
|---|---:|---:|---:|---:|---:|
| Proposed H-DAB dominance | 1.483 | 0.097 | 0.000 | 0.000 | 1.973 |
| CCE-style auto seeds | 1.403 | 0.118 | 0.000 | 0.000 | 1.933 |
| Fixed H-DAB deconvolution | 1.159 | 0.293 | 0.000 | 0.000 | 1.814 |
| OD k-means | 0.335 | 0.095 | 0.000 | 0.000 | 0.549 |
| RGB k-means | 0.314 | 0.120 | 0.000 | 0.000 | 0.546 |

结论：

```text
在三张真实图的无标注 QC 指标上，本插件优于测试的传统算法。
```

## 10. 当前限制

### 10.1 没有人工 ground truth

三张真实图没有人工标注 DAB-positive 区域，因此不能计算真实 F1、IoU、accuracy。

当前真实图比较使用的是无标注 QC：

```text
白背景误判
黑色杂质误判
DAB-positive 与 non-positive 的 DAB 强度分离
Hematoxylin 泄漏
Residual 泄漏
```

### 10.2 当前插件偏高可信阳性

当前策略使用 DAB-dominant 第 85 百分位作为阈值下限。优点是背景干净、淡背景不容易误判；缺点是弱阳性可能被漏掉。

如果后续任务要检测弱阳性，应增加模式选项：

```text
Strict mode: 当前版本，高可信阳性
Sensitive mode: 降低 DAB cloud floor，提高弱阳性召回
Manual mode: 用户手动设置 DAB threshold / dominance threshold
```

### 10.3 Stain vector 固定

当前使用固定 H-DAB vectors。如果不同批次染色颜色偏差很大，建议增加 ROI 校准：

```text
H-only ROI
DAB-strong ROI
background ROI
```

## 11. 后续工程建议

建议下一版增加：

```text
1. Strict / Balanced / Sensitive 三档模式
2. ROI-based stain vector calibration
3. Batch folder processing UI
4. Overlay preview image
5. 可选择只输出三张 color deconvolution 图
6. 保存参数 JSON
7. 与 Fiji Results Table 更深集成
```

## 12. 验收标准

当前版本满足以下验收项：

```text
插件已安装到 /Users/tom/Fiji/plugins
Fiji 菜单可识别 plugins.config
可处理 RGB H-DAB 图像
输出三张彩色 color deconvolution 图
输出三张浓度图
输出 DAB-positive mask
输出 measurement CSV
已在三张真实图上运行
已与传统方法做无标注 QC 比较
真实图 QC score 高于测试的传统方法
```

未满足或需要外部输入的项：

```text
真实图 formal F1 / IoU / accuracy
```

原因：

```text
缺少人工标注或病理专家 ROI。
```
