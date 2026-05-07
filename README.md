# toy-LCAO-bandstructure

+ 这是一个极简模型，通过紧束缚方法/原子轨道线性组合法机制来直接计算能带结构。

+ 脚本在 Matlab 上运行。

+ 需要自行设定原子轨道和原子势场的形式。示例中给的是简单的单个斯莱特轨道+库仑屏蔽势。

+ 可以绘制 k 平面上的能带图，也可以绘制沿指定高对称路径的能带图（需要自行输入节点坐标）。

+ 该脚本**并非**为科研或工业计算而设计，编写它主要是为实践我在本科固体物理课程中学到的紧束缚方法。但它确实可以计算很多种结构。我在我的本科毕业设计中用它复现了Si-mp149, GaAs-mp2534, C-mp48的能带。还可以探索更多。

This is a minimal model based on TBM/LCAO mechanism used for electron bands calculation. The model features direct integral computaion and runs on Matlab. Atomic orbitals and single atomic potentials are required. In examples we assigned single Slater-type orbitals and Coulomb screening potentials. 

The scripts are NOT designed for research or industrial-level calculations. It is a toy model that I developed to test the tight-binding method I learned during my undergraduate solid states physics course. But the LCAO mechanism is clear and stressed (which makes it suitable for bands displaying and other educational purposes). It is also involved in my graduation project.

So far the model has replicated bands patterns of Si-mp149, GaAs-mp2534 and C-mp48. There are more to be discovered.
