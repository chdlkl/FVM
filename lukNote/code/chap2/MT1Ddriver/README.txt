在IDE下运行此程序
1. win32做如下设置：
打开工程属性->Fortran->Labraries->Use Intel Math Kernel Library设置为Sequential
打开工程属性->Linker->Input->Additional Dependencies中写入mkl_lapack95.lib

2. x64做如下设置
打开工程属性->Fortran->Labraries->Use Intel Math Kernel Library设置为Sequential
打开工程属性->Linker->Input->Additional Dependencies中写入mkl_lapack95_lp64.lib

如果遇到stack overflow运行时错误，做如下设置
打开工程属性->Fortran->Linker->System中的Stack Reserve Size与Stack Commit Size均设置为1000000000
打开工程属性->Fortran->Linker->System中的Enable Large Addresses设置为Support Addresses Larger Than 2 GB (/LARGEADDRESSAWARE)


命令行运行此程序
1. dos下运行
win32
ifort /Qmkl mkl_lapack95.lib MT1Ddriver.f90 ddx.f90 sdiag.f90 ave.f90 /link /stack:100000000
test1Dem.exe

x64环境
ifort /Qmkl mkl_lapack95_lp64.lib MT1Ddriver.f90 ddx.f90 sdiag.f90 ave.f90 /link /stack:100000000
test1Dem.exe
在代码没有bug的情况下，建议用x64的DOS运行，或者使用x64的release模式运行

linux下与dos下类似，要注意的是，win的lib文件不能用于linux，linux的Intel fortran的函数库有自带的lib文件
