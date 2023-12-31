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
![2Dplot](https://github.com/Yangjianzin/EZDrawing/assets/22924622/76b144d7-8c8e-4122-bcd7-cdb7670b16f6)

## 如何開始繪製一個3D plot
```C sharp
EZDrawing Viewer = new EZDrawing(5, 200, ViewerType.Draw3D);//建立一個3D畫布

Viewer.backgroundColor = Color.White;

Viewer.InitViewer();

for (double deg = 0; deg <= 360; deg += 0.5)//在高度0mm,0.5mm,1mm各畫一個半徑1mm的紅色圓
{
    D2Point p = D2Math.getCirclePoint(1, deg);
    Viewer.DrawPoint(new Point3D(p.X, p.Y,0), 0.05, Color.Red);
    Viewer.DrawPoint(new Point3D(p.X, p.Y, 0.5), 0.05, Color.Red);
    Viewer.DrawPoint(new Point3D(p.X, p.Y, 1), 0.05, Color.Red);
}

ViewerImage.Image = Viewer.ViewerImage;
```
![3DPlot](https://github.com/Yangjianzin/EZDrawing/assets/22924622/786b14a5-d304-472d-8535-519089a253e6)
