========================================================================
     Water flow driven Salient Object Detection at 180 fps
========================================================================
* Copyright(c) 2016-2017 Xiaoming Huang, <huangxm14@mails.tsinghua.edu.cn>.
* Last update: 17-01-2017. 

This is the demo of our work "Water flow driven Salient Object Detection at 180 fps",
Note that our demo program developed with visual studio 2012 and OpenCV 2.4.9, we have tested the program on a 64-bit PC with Win7 OS,
we achieve 180 fps speed performance on Intel Core i5-4590 CPU @ 3.3 GHz and 8GB RAM.
If you have any problem, please don't hesitate to contact with huangxm14@mails.tsinghua.edu.cn.

Note that we find that time cost of demo program fluctuate due to some reason, we suggest run 5 times and take average time cost.


1. Directory introduction
--ECCSD: 100 images and ground-truth selected from ECCSD dataset, saliency result of ours, MST [45] and MB+ [40];
--MSRA: 100 images and ground-truth selected from MSRA10K dataset, saliency result of ours, MST [45] and MB+ [40];
--Demo.bat: one .bat program to detect saliency of ECCSD and MSRA directory
--SaliencyWater.exe: our program

2. Run Demo Program
You can run Demo.bat, input and output files are in directory Dataset.
