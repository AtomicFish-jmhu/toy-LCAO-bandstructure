# toy-LCAO-bandstructure

+ 这是一个极简模型，通过紧束缚方法/原子轨道线性组合法机制来直接计算能带结构。

+ 脚本在 Matlab 上运行。

+ 需要自行设定原子轨道和原子势场的形式。示例中给的是简单的单个斯莱特轨道+库仑屏蔽势。目前程序只接收球对称标势，不支持非局域势。

+ 可以绘制 k 平面上的能带图，也可以绘制沿指定高对称路径的能带图（需要自行输入节点坐标）。

+ 编写它主要是为实践我在本科固体物理课程中学到的紧束缚方法，适于演示目的。但它确实可以计算很多种结构，尤其是几何机制占主导而对势场参数不敏感的情况下。我在我的本科毕业设计中用它复现了Si-mp149, GaAs-mp2534, C-mp48等物质的能带。


This is a minimal model based on TBM/LCAO mechanism used for bands structure calculation. The model features direct integral computaion and runs on Matlab. The model inputs self-assigned atomic orbitals and single atomic potentials. Single Slater-type orbitals and Coulomb screening potentials are adopted in the examples. Currently nonlocal potentials are not supported.

The scripts is a toy model that I developed to test the tight-binding method I learned during my undergraduate solid states physics course. But the LCAO mechanism is clear and stressed (which makes it suitable for bands demonstration and other education purposes). The model is suitable for materials where the topology dominates and that are not sensitive to parameteres. It is also involved in my graduation project.

So far the model has replicated band patterns of Si-mp149, GaAs-mp2534 and C-mp48.
