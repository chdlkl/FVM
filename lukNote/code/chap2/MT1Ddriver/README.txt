��IDE�����д˳���
1. win32���������ã�
�򿪹�������->Fortran->Labraries->Use Intel Math Kernel Library����ΪSequential
�򿪹�������->Linker->Input->Additional Dependencies��д��mkl_lapack95.lib

2. x64����������
�򿪹�������->Fortran->Labraries->Use Intel Math Kernel Library����ΪSequential
�򿪹�������->Linker->Input->Additional Dependencies��д��mkl_lapack95_lp64.lib

�������stack overflow����ʱ��������������
�򿪹�������->Fortran->Linker->System�е�Stack Reserve Size��Stack Commit Size������Ϊ1000000000
�򿪹�������->Fortran->Linker->System�е�Enable Large Addresses����ΪSupport Addresses Larger Than 2 GB (/LARGEADDRESSAWARE)


���������д˳���
1. dos������
win32
ifort /Qmkl mkl_lapack95.lib MT1Ddriver.f90 ddx.f90 sdiag.f90 ave.f90 /link /stack:100000000
test1Dem.exe

x64����
ifort /Qmkl mkl_lapack95_lp64.lib MT1Ddriver.f90 ddx.f90 sdiag.f90 ave.f90 /link /stack:100000000
test1Dem.exe
�ڴ���û��bug������£�������x64��DOS���У�����ʹ��x64��releaseģʽ����

linux����dos�����ƣ�Ҫע����ǣ�win��lib�ļ���������linux��linux��Intel fortran�ĺ��������Դ���lib�ļ�
