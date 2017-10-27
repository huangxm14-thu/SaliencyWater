================================================================================================
               Xiaoming Huang, and Yu-Jin Zhang,
               Water flow driven salient object detection at 180 fps,
               Pattern Recognition
================================================================================================
* Copyright(c)2015-2017 Xiaoming Huang, <huangxm14@mails.tsinghua.edu.cn>.
* Last update:26-10-2017. 

This is the code of above work. The source code is developed with visual studio 2012, intel ipp library and opencv 2.4.9.
We have test the program on a 64-bit PC with Win7 OS,achieve 180 FPS speed performance on Intel Core i5-4590 CPU @ 3.3 GHz and 8GB RAM.
If you have any problem, please don't hesitate to contact with huangxm14@mails.tsinghua.edu.cn.

***Note that we find that time cost of program fluctuate due to some reason, we suggest run 5 times and take average time cost.


1.Directory introduction
--Demo: One demo program
--SaliencyWater.sln: solution file (open with visual studio 2012)
--opencv2.4.9:header files and lib of opencv
--test.cpp: test file of our method
--Saliency.cpp: detail implementation of our method

2.Build Project
Open solution file saliencyWater.sln with VS2012, 
modify header file path of intel ipp library in saliency.cpp
select Release @ x64 mode, build.

3.Run
You need set directory of input images, e.g., .\Demo\MSRA, then run the program, the final result will be also saved in the directory.



