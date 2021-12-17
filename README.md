# LiverVesselSegmentMRF
肝脏血管的去噪，特征提取，MRF分割方法
liver vessel denois, feature map extract and segment by MRF method


## 使用方法(Usage)
解压代码文件夹中的libraryyoumayuses .7z文件。如果你在自己的论文中使用它们，你应该引用它们的作者。  
更改文件main.m中标记的路径，运行main.m。代码将自动计算一些度量，您可以在工作空间中找到，以单词“crit_”开头。  
在文件夹“result/”中检查结果。特征图导出到“feature_map/”文件夹  
Unzip the LibraryYouMayUsed.7z file in the code folder. You should cite their author if you use them in you own paper.  
Change the path marked in the main.m file, and run. The code will automatically calculate some metric, which you can find in workspace start with the words "crit_".  
Check the result in the folder names "result/". Feature map may be exported in the folder names "feature_map/".

## DIY your own dataset
You can use your own dataset by simply changing the function start with the words "DataLoader" just like I did. 
Though you have to provide mask file of liver, regrettably. Label is not exactly needed, but that may cause some mistake and you have to debug.  
