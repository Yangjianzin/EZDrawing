# EZDrawing
Use C# to Draw 2D and 3D plot
## 在使用前請先將EZDrawing加入專案，並加入參考
## 在程式碼上方加入
```C sharp
using EZDrawingDLL;
```
## 如何開始繪製一個2D plot
```C sharp
EZDrawing Viewer = new EZDrawing(5, 200, ViewerType.CenterXY);//建立一個bitmap 1200x1200 畫布

Viewer.backgroundColor = Color.White;//設定背景顏色

Viewer.InitViewer();//先初始化畫布

Viewer.DrawCenterCross(1, Color.Black);//繪製一個XY坐標系

Viewer.DrawPoint(1, 1, 0.1, Color.Red);//在位置X=1mm,Y=1mm畫一個半徑為0.1mm的紅點

ViewerImage.Image = Viewer.ViewerImage;//將bitmap顯示在picturebox上
```
