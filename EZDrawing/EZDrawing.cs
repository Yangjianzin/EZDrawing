using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Runtime.InteropServices;
using System.Drawing.Imaging;


namespace EZDrawingDLL
{
    public enum ViewerType
    {
        CenterXY, X_ZDepth, Y_XTime, Draw3D, Infomation
    }
    public class D2Point : IDisposable
    {
        public double X;

        public double Y;
        
        public int Index;

        public double Z;
        
        public bool IsBufferPoint = false;

        public double Angle = 0;
        public void Dispose()
        {
        }

        public D2Point()
        {

        }

        public D2Point(double x, double y) : this(x, y, 0) { }

        public D2Point(double x, double y, double z, bool buffer = false)
        {
            X = x;
            Y = y;
            Index = 5;//預設是開
            Z = z;
            IsBufferPoint = buffer;
        }

        public D2Point Clone()
        {
            return (D2Point)this.MemberwiseClone();
        }

        public static bool TryParse(string str, out D2Point p)
        {
            p = new D2Point();
            string[] data = str.Split(',');
            if (data.Length == 4)
            {
                if (!double.TryParse(data[0], out p.X))
                    return false;
                if (!double.TryParse(data[1], out p.Y))
                    return false;
                if (!int.TryParse(data[2], out p.Index))
                    return false;
                if (!double.TryParse(data[3], out p.Z))
                    return false;
                if (p.Index == 5)
                    p.IsBufferPoint = true;
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    public class D2Math
    {
        public const double M_PI = 3.1415926535;
        /// <summary>
        /// 計算允許最小誤差值
        /// </summary>
        public const double AllowErrorValue = 0.0001;
        /// <summary>
        /// Deg 轉換 Rad
        /// </summary>
        /// <param name="degrees"></param>
        /// <returns></returns>
        public static double ConvertDegreesToRadians(double degrees)
        {
            double radians = (double)(M_PI / 180.0) * (double)degrees;
            return (radians);
        }
        /// <summary>
        /// 取得圓上的點
        /// </summary>
        /// <param name="r"></param>
        /// <param name="degrees"></param>
        /// <returns></returns>
        public static D2Point getCirclePoint(double r, double degrees)
        {
            double dx = ((double)Math.Cos((double)ConvertDegreesToRadians(degrees)) * r);
            double dy = (-(double)Math.Sin((double)ConvertDegreesToRadians(degrees)) * r);
            dx = Math.Round(dx, 8);
            dy = Math.Round(dy, 8);
            return new D2Point(dx, dy);
        }
        /// <summary>
        /// 取得圓上的點
        /// </summary>
        /// <param name="r"></param>
        /// <param name="degrees"></param>
        /// <returns></returns>
        public static D2Point getOvalPoint(double r1, double r2, double degrees)
        {
            double dx = ((double)Math.Cos((double)ConvertDegreesToRadians(degrees)) * r1);
            double dy = (-(double)Math.Sin((double)ConvertDegreesToRadians(degrees)) * r2);
            return new D2Point(dx, dy);
        }
        /// <summary>
        /// 計算距離
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double getDistance(D2Point p1, D2Point p2)
        {
            double dx = p1.X - p2.X;
            double dy = p1.Y - p2.Y;
            double line = (double)Math.Sqrt((double)(dx * dx + dy * dy));
            return line;
        }
        public static double getDistance1(D2Point p1, D2Point p2)
        {
            double dx = p1.X - p2.X;
            double dy = p1.Y - p2.Y;
            return (dx * dx + dy * dy);
          
        }

        /// <summary>  
        /// 判断点是否在多边形内.  
        /// ----------原理----------  
        /// 注意到如果从P作水平向左的射线的话，如果P在多边形内部，那么这条射线与多边形的交点必为奇数，  
        /// 如果P在多边形外部，则交点个数必为偶数(0也在内)。  
        /// </summary>  
        /// <param name="checkPoint">要判断的点</param>  
        /// <param name="polygonPoints">多边形的顶点</param>  
        /// <returns></returns>  
        public static bool IsInPolygon(D2Point checkPoint, List<D2Point> polygonPoints)
        {

            bool inside = false;

            int pointCount = polygonPoints.Count;
            D2Point p1, p2;
            for (int i = 0, j = pointCount - 1; i < pointCount; j = i, i++)//第一个点和最后一个点作为第一条线，之后是第一个点和第二个点作为第二条线，之后是第二个点与第三个点，第三个点与第四个点...  
            {
                p1 = polygonPoints[i];
                p2 = polygonPoints[j];
                if (checkPoint.Y < p2.Y)
                {//p2在射线之上  
                    if (p1.Y <= checkPoint.Y)
                    {//p1正好在射线中或者射线下方  
                        if ((checkPoint.Y - p1.Y) * (p2.X - p1.X) > (checkPoint.X - p1.X) * (p2.Y - p1.Y))//斜率判断,在P1和P2之间且在P1P2右侧  
                        {
                            //射线与多边形交点为奇数时则在多边形之内，若为偶数个交点时则在多边形之外。  
                            //由于inside初始值为false，即交点数为零。所以当有第一个交点时，则必为奇数，则在内部，此时为inside=(!inside)  
                            //所以当有第二个交点时，则必为偶数，则在外部，此时为inside=(!inside)  
                            inside = (!inside);
                        }
                    }
                }
                else if (checkPoint.Y < p1.Y)
                {
                    //p2正好在射线中或者在射线下方，p1在射线上  
                    if ((checkPoint.Y - p1.Y) * (p2.X - p1.X) < (checkPoint.X - p1.X) * (p2.Y - p1.Y))//斜率判断,在P1和P2之间且在P1P2右侧  
                    {
                        inside = (!inside);
                    }
                }
            }
            return inside;
        }

        public static double Sqrt(double value)
        {
            return (double)(Math.Sqrt((double)value));
        }

        public static double Sin(double value)
        {
            return (double)(Math.Sin((double)value));
        }

        public static double Cos(double value)
        {
            return (double)(Math.Cos((double)value));
        }
        public static double Acos(double value)
        {
            return (double)(Math.Acos((double)value));
        }
        public static double ConvertRadiansToDegrees(double Radians)
        {
            double Degree = Radians * (double)(180.0 / M_PI);
            return Degree;
        }
    }
    public class Point3D
    {
        public double X;

        public double Y;

        public double Z;

        public Point3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }

        public Point3D()
        {
            // TODO: Complete member initialization
        }

        static Point3D EulerRotationZ(double x, double y, double z, int degree)
        {
            Point3D p = new Point3D();
            p.X = x * D2Math.Cos(D2Math.ConvertDegreesToRadians(degree)) - y * Math.Sin(D2Math.ConvertDegreesToRadians(degree));
            p.Y = y * Math.Cos(D2Math.ConvertDegreesToRadians(degree)) + x * Math.Sin(D2Math.ConvertDegreesToRadians(degree));
            p.Z = z;
            return p;
        }

        static Point3D EulerRotationX(double x, double y, double z, int degree)
        {
            Point3D p = new Point3D();
            p.X = x;
            p.Y = y * D2Math.Cos(D2Math.ConvertDegreesToRadians(degree)) - z * D2Math.Sin(D2Math.ConvertDegreesToRadians(degree));
            p.Z = y * D2Math.Sin(D2Math.ConvertDegreesToRadians(degree)) + z * D2Math.Cos(D2Math.ConvertDegreesToRadians(degree));
            return p;
        }

        static Point3D EulerRotationY(double x, double y, double z, int degree)
        {
            Point3D p = new Point3D();
            p.X = x * D2Math.Cos(D2Math.ConvertDegreesToRadians(degree)) + z * D2Math.Sin(D2Math.ConvertDegreesToRadians(degree));
            p.Y = y;
            p.Z = -x * D2Math.Sin(D2Math.ConvertDegreesToRadians(degree)) + z * D2Math.Cos(D2Math.ConvertDegreesToRadians(degree));
            return p;
        }

        public static Point3D TranEulerRotation(Point3D p, EulerRotation ER)
        {

            p = EulerRotationX(p.X, p.Y, p.Z, ER.X);
            p = EulerRotationY(p.X, p.Y, p.Z, ER.Y);
            p = EulerRotationZ(p.X, p.Y, p.Z, ER.Z);
            return p;
        }

    }

    public class EulerRotation
    {
        public int X;

        public int Y;

        public int Z;

        public EulerRotation()
        {

        }

        public EulerRotation(int R)
            : this(R, R, R)
        {

        }

        public EulerRotation(int x, int y, int z)
        {
            X = x;
            Y = y;
            Z = z;
        }


    }

    public class DirectBitmap : IDisposable
    {
        public Bitmap Bitmap { get; private set; }
        public Int32[] Bits { get; private set; }
        public bool Disposed { get; private set; }
        public int Height { get; private set; }
        public int Width { get; private set; }

        protected GCHandle BitsHandle { get; private set; }

        public DirectBitmap(int width, int height)
        {
            Width = width;
            Height = height;
            Bits = new Int32[width * height];
            BitsHandle = GCHandle.Alloc(Bits, GCHandleType.Pinned);
            Bitmap = new Bitmap(width, height, width * 4, PixelFormat.Format32bppPArgb, BitsHandle.AddrOfPinnedObject());
        }

        public void SetPixel(int x, int y, Color colour)
        {
            int index = x + (y * Width);
            int col = colour.ToArgb();

            Bits[index] = col;
        }

        public Color GetPixel(int x, int y)
        {
            int index = x + (y * Width);
            int col = Bits[index];
            Color result = Color.FromArgb(col);

            return result;
        }

        public void Dispose()
        {
            if (Disposed) return;
            Disposed = true;
            Bitmap.Dispose();
            BitsHandle.Free();
        }
    }

    public class EZDrawing:IDisposable
    {
        public DirectBitmap _Image;
        public Bitmap ViewerImage
        {
            get
            {
                return _Image.Bitmap;
            }
         
        }
        //中心點座標
        public PointF Center;

        private Graphics g;
        //畫布背景顏色
        public Color backgroundColor = Color.Black;

        public int ViewerSize = 1000;
        //1 unit spilt ?
        public int Resolution = 100;

        private int InfoAreaSize_mm = 1;

        private int DrawAreaSize_mm = 5;

        public double MaxSize
        {
            get
            {
                return ViewerSize / Resolution / 2.0;
            }
        }

        private double _FontSize = 0;

        public bool IsFasterDraw = false;
        public double FontSize
        {
            get
            {
                return _FontSize;
            }
        }
        
        private Font defalutFont = new Font("Tahoma", 25F, FontStyle.Regular, GraphicsUnit.Point, ((byte)(0)));

        private ViewerType Vtype = ViewerType.CenterXY;

        public EZDrawing()
        {

            Center = new PointF(ViewerSize / 2, ViewerSize / 2);
            _Image = new DirectBitmap(ViewerSize, ViewerSize);
            AutoAdjustFont(0.5);
            InitViewer();
        }

        public EZDrawing(int DrawAreaSize_mm, int resolution_In_1mm, ViewerType vt = ViewerType.CenterXY, double FontSize = 0.3)
        {
            this.DrawAreaSize_mm = DrawAreaSize_mm;
            ViewerSize = (DrawAreaSize_mm + InfoAreaSize_mm) * 2 * resolution_In_1mm;
            Resolution = resolution_In_1mm;
            Center = new PointF(ViewerSize / 2, ViewerSize / 2);
            _Image = new DirectBitmap(ViewerSize, ViewerSize);
            Vtype = vt;
            AutoAdjustFont(FontSize);
            InitViewer();
            Zbuffer = new float[ViewerSize, ViewerSize];
            for (int i = 0; i < ViewerSize; i++)
                for (int j = 0; j < ViewerSize; j++)
                    Zbuffer[i, j] = float.MinValue;
        }

        public void InitViewer()
        {
            try
            {
                g = Graphics.FromImage(ViewerImage);
                g.Clear(backgroundColor);
                g.SmoothingMode = SmoothingMode.AntiAlias;

                g.InterpolationMode = InterpolationMode.HighQualityBilinear; //設定高品質插值法
                g.SmoothingMode = SmoothingMode.HighQuality; //設定高品質,低速度呈現平滑程度
                g.CompositingQuality = CompositingQuality.HighQuality;
                DrawAreaOutLine(DrawAreaSize_mm, 0.01, Color.Gray);
            }
            catch
            {

            }
        }

        public Color Rainbow(float progress)
        {

            progress /= 1.1f;
            float div = (Math.Abs(progress % 1) * 6);
            int ascending = (int)((div % 1) * 255);
            int descending = 255 - ascending;

            switch ((int)div)
            {
                case 0:
                    return Color.FromArgb(255, 255, ascending, 0);
                case 1:
                    return Color.FromArgb(255, descending, 255, 0);
                case 2:
                    return Color.FromArgb(255, 0, 255, ascending);
                case 3:
                    return Color.FromArgb(255, 0, descending, 255);
                case 4:
                    return Color.FromArgb(255, ascending, 0, 255);
                default://case 5:
                    return Color.FromArgb(255, 255, 0, descending);
            }
        }

        public void DrawColorTable(double MaxValue, int ColorSpace = 5)
        {
            int MaxIndex = (int)((MaxValue % ColorSpace == 0) ? MaxValue / ColorSpace : (MaxValue + ColorSpace) / ColorSpace);
            double P = DrawAreaSize_mm + InfoAreaSize_mm;
            for (int c = 0; c <= MaxIndex; c++)
            {
                float v = c * (float)MaxValue / (float)MaxIndex;
                DrawString(v.ToString("0.00"), -P + 0.5, DrawAreaSize_mm - 0.2 - c / 3f, null, Rainbow(v / (float)MaxValue));
            }
        }

        public double offsetY
        {
            get
            {
                double offsetY = 0;
                if (Vtype == ViewerType.X_ZDepth)
                    offsetY = -DrawAreaSize_mm;

                return offsetY;
            }
        }

        public double offsetX
        {
            get
            {
                double offsetX = 0;
                if (Vtype == ViewerType.Y_XTime)
                    offsetX = -DrawAreaSize_mm;
                return offsetX;
            }
        }

        //Float RainBowColor
        public void DrawPoint(double x, double y, double r, float Rainbowc)
        {
            DrawPoint(x, y, r, Rainbow(Rainbowc));
        }

        public void DrawPoint(Point3D p, double r, float Rainbowc)
        {
            DrawPoint(p, r, Rainbow(Rainbowc));
        }

        public void DrawCircle(double x, double y, double r, double w, float Rainbowc)
        {
            DrawCircle(x, y, r, w, Rainbow(Rainbowc));
        }

        public void DrawPoint(Point3D p, EulerRotation ER, double r, float Rainbowc)
        {
            DrawPoint(p, ER, r, Rainbow(Rainbowc));
        }

        public void DrawArc(double x, double y, double r, double startAngle, double endAngle, double w, float Rainbowc)
        {
            DrawArc(x, y, r, startAngle, endAngle, w, Rainbow(Rainbowc));
        }

        public void DrawString(string str, double x, double y, Font f, float Rainbowc)
        {
            DrawString(str, x, y, f, Rainbow(Rainbowc));
        }


        //Color
        public void DrawPoint(double x, double y, double r, Color color)
        {
            DrawBasePoint((float)(x + offsetX), (float)(y + offsetY), (float)r, color);
        }

        public void DrawPoint(Point3D p, double r, Color color)
        {
            switch (Vtype)
            {
                case ViewerType.CenterXY:
                    DrawBasePoint((float)(p.X + offsetX), (float)(p.Y + offsetY), (float)r, color);
                    break;

                case ViewerType.X_ZDepth:
                    DrawBasePoint((float)(p.X + offsetX), (float)(p.Z + offsetY), (float)r, color);
                    break;

                case ViewerType.Draw3D:
                    DrawPoint(p, new EulerRotation(320), r, color);
                    break;
            }

        }

        public void DrawPoint(Point3D p, EulerRotation ER, double r, Color color)
        {
            Point3D p3d = Point3D.TranEulerRotation(p, ER);
            double offsetY3D = -DrawAreaSize_mm * 0.5;
            DrawBase3DPoint((float)(p3d.X + offsetX), (float)(p3d.Y + offsetY + offsetY3D), (float)p3d.Z, (float)r, color);
        }

        public void DrawCircle(double x, double y, double r, double w, Color color)
        {
            DrawBaseCircle((float)(x + offsetX), (float)(y + offsetY), (float)r, (float)w, color);
        }
        public void DrawEllipse(double x, double y, double r, double r1, double w, Color color)
        {
            DrawBaseEllipse((float)(x + offsetX), (float)(y + offsetY), (float)r, (float)r1, (float)w, color);
        }
        public void DrawArc(double x, double y, double r, double startAngle, double endAngle, double w, Color color, bool IsDashLine = false)
        {
            DrawBaseArc((float)(x + offsetX), (float)(y + offsetY), (float)r, (float)startAngle, (float)endAngle, (float)w, color, IsDashLine);
        }

        public void DrawString(string str, double x, double y, Font f, Color color, bool IsLeft = false)
        {
            if (f == null)
            {
                f = defalutFont;
            }
            DrawBaseString(str, (float)(x + offsetX), (float)(y + offsetY), f, color, IsLeft);
        }

        public void DrawInfomation(string str, int lines, Color color, bool IsLeft = false)
        {
            double dx = -(DrawAreaSize_mm + InfoAreaSize_mm);
            double dy = (DrawAreaSize_mm + InfoAreaSize_mm - (lines + 1) * _FontSize);
            if (Vtype == ViewerType.X_ZDepth)
            {
                dy = ((DrawAreaSize_mm * 2 + InfoAreaSize_mm) - (lines + 1) * _FontSize);
            }
            DrawString(str, dx, dy, null, color, IsLeft);

        }

        public void DrawLine(double x1, double y1, double x2, double y2, double width, Color color, bool IsDashLine = false)
        {
            DrawBaseLine((float)(x1 + offsetX), (float)(y1 + offsetY), (float)(x2 + offsetX), (float)(y2 + offsetY), (float)width, color, IsDashLine);
        }

        public void DrawLine(double x1, double y1, double x2, double y2, double width, float Rainbowc, bool IsDashLine = false)
        {
            DrawLine(x1, y1, x2, y2, width, Rainbow(Rainbowc), IsDashLine);
        }


        //外框方法
        public void DrawAreaOutLine(double AreaLength, double outline, Color color)
        {
            ViewerType temp = Vtype;
            Vtype = ViewerType.CenterXY;
            DrawLine(AreaLength, AreaLength, -AreaLength, AreaLength, outline, color, true);
            DrawLine(-AreaLength, -AreaLength, AreaLength, -AreaLength, outline, color, true);
            DrawLine(AreaLength, -AreaLength, AreaLength, AreaLength, outline, color, true);
            DrawLine(-AreaLength, -AreaLength, -AreaLength, AreaLength, outline, color, true);
            Vtype = temp;
        }

        public void Draw_3D_Box(EulerRotation ER)
        {
            double dx = -DrawAreaSize_mm / 2;
            double dy = -DrawAreaSize_mm / 2;
            double dr = DrawAreaSize_mm;

            Point3D p0 = new Point3D(dx, dy, 0);
            Point3D p1 = new Point3D(dx, dy + dr, 0);
            Point3D p2 = new Point3D(dx + dr, dy, 0);
            Point3D p3 = new Point3D(dx, dy, dr);

            Point3D p00 = Point3D.TranEulerRotation(p0, ER);
            Point3D p10 = Point3D.TranEulerRotation(p1, ER);
            Point3D p20 = Point3D.TranEulerRotation(p2, ER);
            Point3D p30 = Point3D.TranEulerRotation(p3, ER);

            DrawString("x", p10.X, p10.Y, null, Color.Red, true);
            DrawString("y", p20.X, p20.Y, null, Color.Red, true);
            DrawString("z", p30.X, p30.Y, null, Color.Red, true);

            DrawLine(p00.X, p00.Y, p10.X, p10.Y, 0.02, Color.Gray);
            DrawLine(p00.X, p00.Y, p20.X, p20.Y, 0.02, Color.Gray);
            DrawLine(p00.X, p00.Y, p30.X, p30.Y, 0.02, Color.Gray);
        }

        public void Draw_X_Z_Box(double ZSpace, double Zfactor, double XSpace, Color color, Font font = null, bool IsShowUnit = true)
        {

            double LineWidth = 0.015;
            double uLength = 0.2;
            if (font == null)
            {
                font = defalutFont;
            }
            DrawLine(-DrawAreaSize_mm, 0, -DrawAreaSize_mm, DrawAreaSize_mm * 2, LineWidth, color);
            {
                double eachZSpace = DrawAreaSize_mm * 2 / ZSpace;
                for (int i = 0; i < eachZSpace; i++)
                {
                    double tSpace = i * ZSpace;
                    double dstrH = _FontSize + uLength;
                    double dstrW = _FontSize + uLength;
                    DrawLine(-DrawAreaSize_mm, tSpace, -DrawAreaSize_mm + uLength, tSpace, LineWidth, color);
                    double minSpace = ZSpace / 10;
                    for (int j = 0; j < 10; j++)
                    {
                        DrawLine(-DrawAreaSize_mm, tSpace + j * minSpace, -DrawAreaSize_mm + uLength / 2, tSpace + j * minSpace, LineWidth, color);
                    }
                    if (IsShowUnit)
                    {

                        string str = string.Format("{0:0}", tSpace / Zfactor);
                        DrawString(str, -DrawAreaSize_mm - dstrH, tSpace, font, color);
                    }
                }
            }
            {
                double eachXSpace = DrawAreaSize_mm / XSpace;
                for (int i = 0; i < eachXSpace; i++)
                {
                    double tSpace = i * XSpace;
                    double dstrH = _FontSize;
                    DrawLine(tSpace, uLength, tSpace, 0, LineWidth, color);
                    DrawLine(-tSpace, uLength, -tSpace, 0, LineWidth, color);
                    if (IsShowUnit)
                    {
                        DrawString((-tSpace).ToString(), -tSpace, -dstrH, font, color);
                        DrawString((tSpace).ToString(), tSpace, -dstrH, font, color);
                    }
                }
            }
        }

        public void Draw_Y_T_Box(double YSpace, double Yfactor, double TimeSpace, double Timefactor, Color color, Font font = null, bool IsShowUnit = true)
        {

            double LineWidth = 0.015;
            double uLength = 0.2;
            if (font == null)
            {
                font = defalutFont;
            }
            //DrawLine(-DrawAreaSize_mm, 0, -DrawAreaSize_mm, DrawAreaSize_mm * 2, LineWidth, color);
            {
                double eachZSpace = DrawAreaSize_mm / YSpace;
                for (int i = 0; i < eachZSpace; i++)
                {
                    double tSpace = i * YSpace;
                    double dstrH = _FontSize + uLength;
                    double dstrW = _FontSize + uLength;
                    DrawLine(0, tSpace, uLength, tSpace, LineWidth, color);
                    DrawLine(0, -tSpace, uLength, -tSpace, LineWidth, color);
                    double minSpace = YSpace / 10;
                    for (int j = 0; j < 10; j++)
                    {
                        DrawLine(0, tSpace + j * minSpace, uLength / 2, tSpace + j * minSpace, LineWidth, color);
                        DrawLine(0, -tSpace - j * minSpace, uLength / 2, -tSpace - j * minSpace, LineWidth, color);

                    }
                    if (IsShowUnit)
                    {
                        string str = string.Format("{0:0}", tSpace / Yfactor);
                        DrawString(str, 0 - dstrH, tSpace, font, color);
                        string str1 = string.Format("{0:0}", -tSpace / Yfactor);
                        DrawString(str1, 0 - dstrH, -tSpace, font, color);
                    }
                }
            }
            {
                double eachXSpace = DrawAreaSize_mm * 2 / TimeSpace;
                for (int i = 0; i <= eachXSpace; i++)
                {
                    double tSpace = i * TimeSpace;
                    double dstrH = _FontSize;

                    DrawLine(tSpace, DrawAreaSize_mm - uLength, tSpace, DrawAreaSize_mm, LineWidth, color);

                    DrawLine(tSpace, -DrawAreaSize_mm + uLength, tSpace, -DrawAreaSize_mm, LineWidth, color);
                    if (IsShowUnit)
                    {
                        string str = string.Format("{0:0}", tSpace / Timefactor);
                        DrawString(str, tSpace, -DrawAreaSize_mm + -dstrH, font, color);
                    }
                }
            }
        }

        public void DrawCenterCross(double UnitSpace, Color color, Font font = null, bool IsShowUnit = true)
        {
            double DrawLength = DrawAreaSize_mm;
            double LineWidth = 0.02;
            double eachSpace = DrawLength / UnitSpace;
            double uLength = 0.1;
            if (font == null)
            {
                font = defalutFont;
            }

            DrawLine(DrawLength, 0, -DrawLength, 0, LineWidth, color);
            DrawLine(0, DrawLength, 0, -DrawLength, LineWidth, color);
            for (int i = 0; i < eachSpace; i++)
            {
                double tSpace = i * UnitSpace;
                double dstrH = _FontSize + uLength;
                double dstrW = _FontSize + uLength;
                DrawLine(-uLength, tSpace, uLength, tSpace, LineWidth, color);
                DrawLine(-uLength, -tSpace, uLength, -tSpace, LineWidth, color);

                DrawLine(tSpace, -uLength, tSpace, uLength, LineWidth, color);
                DrawLine(-tSpace, -uLength, -tSpace, uLength, LineWidth, color);
                if (IsShowUnit)
                {
                    if (i != 0)
                    {
                        DrawString((tSpace).ToString(), dstrH - uLength, tSpace, font, color);
                        DrawString((-tSpace).ToString(), dstrH - uLength, -tSpace, font, color);
                        DrawString((tSpace).ToString(), tSpace, dstrW - uLength, font, color);
                        DrawString((-tSpace).ToString(), -tSpace, dstrW - uLength, font, color);
                    }
                }
            }
        }
        public void DrawCenterLine(double UnitSpace, Color color, Font font = null, bool IsShowUnit = true)
        {
            double DrawLength = DrawAreaSize_mm;
            double LineWidth = 0.02;
            double eachSpace = DrawLength / UnitSpace;
            double uLength = 0.1;
            if (font == null)
            {
                font = defalutFont;
            }

            DrawLine(DrawLength, 0, -DrawLength, 0, LineWidth, color);


            DrawLine(DrawLength, -DrawLength, -DrawLength, -DrawLength, LineWidth, color);
            //DrawLine(DrawLength-0.5, DrawLength, DrawLength-0.5, -DrawLength, LineWidth, color);
            for (int i = 0; i < eachSpace; i++)
            {
                double tSpace = i * UnitSpace;
                double dstrH = _FontSize + uLength;
                double dstrW = _FontSize + uLength;
                //DrawLine(-uLength+ DrawLength - 0.5, tSpace, uLength+ DrawLength - 0.5, tSpace, LineWidth, color);
                //DrawLine(-uLength+ DrawLength - 0.5, -tSpace, uLength+ DrawLength - 0.5, -tSpace, LineWidth, color);

                DrawLine(tSpace, -DrawLength - uLength, tSpace, -DrawLength+ uLength, LineWidth, color);
                DrawLine(-tSpace, -DrawLength - uLength, -tSpace, -DrawLength+uLength, LineWidth, color);


                DrawLine(tSpace, -uLength, tSpace, uLength, LineWidth, color);
                DrawLine(-tSpace, -uLength, -tSpace, uLength, LineWidth, color);
                if (IsShowUnit)
                {
                    if (i != 0)
                    {
                        //DrawString((tSpace).ToString(), dstrH - uLength+ DrawLength - 0.5, tSpace, font, color);
                        //DrawString((-tSpace).ToString(), dstrH - uLength+ DrawLength - 0.5, -tSpace, font, color);
                        DrawString((tSpace).ToString(), tSpace, dstrW - uLength, font, color);
                        DrawString((-tSpace).ToString(), -tSpace, dstrW - uLength, font, color);

                        DrawString((tSpace).ToString(), tSpace, -DrawLength + dstrW - uLength, font, color);
                        DrawString((-tSpace).ToString(), -tSpace, -DrawLength + dstrW - uLength, font, color);
                    }
                    else
                    {
                        DrawString((tSpace).ToString(), tSpace, dstrW - uLength, font, color);

                        DrawString((tSpace).ToString(), tSpace, -DrawLength + dstrW - uLength, font, color);
                    }
                }
            }
        }

        public float GetRealX(float x)
        {
            return x + Center.X;
        }

        public float GetRealY(float y)
        {
            return Center.X - y;
        }

        private float[,] Zbuffer = null;

        public void DrawZbuffer(Rectangle r, double z, Brush br)
        {

            for (int i = 0; i < r.Width; i++)
            {
                for (int j = 0; j < r.Height; j++)
                {
                    int dx = i + r.X;
                    int dy = j + r.Y;
                    //Avoid the Out of the index
                    if ((dx < ViewerSize) && (dy < ViewerSize) && (dx > 0) && (dy > 0))
                    {
                        if (Zbuffer[dx, dy] < z)
                        {
                            Zbuffer[dx, dy] = (float)z;
                            g.FillEllipse(br, dx, dy, 2, 2);

                        }
                    }
                }
            }
        }

        public void DrawFastZbuffer(Rectangle r, double z, Color cr)
        {
            for (int i = 0; i < r.Width; i++)
            {
                for (int j = 0; j < r.Height; j++)
                {
                    int dx = i + r.X;
                    int dy = j + r.Y;
                    //Avoid the Out of the index
                    if ((dx < ViewerSize) && (dy < ViewerSize) && (dx > 0) && (dy > 0))
                    {
                        if (Zbuffer[dx, dy] < z)
                        {
                            Zbuffer[dx, dy] = (float)z;
                            _Image.SetPixel(dx, dy, cr);
                        }
                    }
                }
            }
        }

        //底層繪圖
        private void DrawBase3DPoint(float x, float y, float z, float r, Color color)
        {
            Brush br = new SolidBrush(color);
            float dx = x * Resolution;
            float dy = y * Resolution;
            float dz = z * Resolution;
            float dr = r * Resolution;

            Rectangle rr = new Rectangle((int)(GetRealX(dx - dr)), (int)(GetRealY(dy + dr)), (int)(dr * 2 < 1 ? 1 : dr * 2), (int)(dr * 2 < 1 ? 1 : dr * 2));
            DrawFastZbuffer(rr, z, color);
            //DrawZbuffer(rr, dz, br);
        }

        private void DrawBasePoint(float x, float y, float r, Color color)
        {

            float dx = x * Resolution;
            float dy = y * Resolution;
            float dr = r * Resolution;
            if (IsFasterDraw)
            {
                _Image.SetPixel((int)GetRealX(dx), (int)GetRealY(dy), color);
            }
            else
            {
                Brush br = new SolidBrush(color);
                RectangleF rf = new RectangleF(GetRealX(dx - dr), GetRealY(dy + dr), dr * 2, dr * 2);
                g.FillEllipse(br, rf);
            }
        }

        private void DrawBaseArc(float x, float y, float Radius, float startAngle, float endAngle, float w, Color color, bool IsDashLine = false)
        {
            float dw = w * Resolution;
            float dx = x * Resolution;
            float dy = y * Resolution;
            float dr = Radius * Resolution;
            Brush br = new SolidBrush(color);
            Pen pen = new Pen(br, dw);
            if (IsDashLine)
                pen.DashStyle = DashStyle.DashDot;
            RectangleF rf = new RectangleF(GetRealX(dx - dr), GetRealY(dy + dr), dr * 2, dr * 2);
            g.DrawArc(pen, rf, startAngle, endAngle - startAngle);
        }

        private void DrawBaseCircle(float x, float y, float Radius, float w, Color color)
        {
            float dw = w * Resolution;
            float dx = x * Resolution;
            float dy = y * Resolution;
            float dr = Radius * Resolution;
            Brush br = new SolidBrush(color);
            Pen pen = new Pen(br, dw);
            //pen.DashStyle = DashStyle.Dash;
            RectangleF rf = new RectangleF(GetRealX(dx - dr), GetRealY(dy + dr), dr * 2, dr * 2);
            g.DrawEllipse(pen, rf);
        }
        private void DrawBaseEllipse(float x, float y, float Radius, float Radius1, float w, Color color)
        {
            float dw = w * Resolution;
            float dx = x * Resolution;
            float dy = y * Resolution;
            float dr = Radius * Resolution;
            float dr1 = Radius1 * Resolution;
            Brush br = new SolidBrush(color);
            Pen pen = new Pen(br, dw);
            //pen.DashStyle = DashStyle.Dash;
            RectangleF rf = new RectangleF(GetRealX(dx - dr), GetRealY(dy + dr1), dr * 2, dr1 * 2);
            g.DrawEllipse(pen, rf);
        }
        private Size getFontSize(string s, Font f)
        {
            Size sss = new Size();
            using (Graphics gr = Graphics.FromImage(new Bitmap(1, 1)))
            {
                SizeF size = gr.MeasureString(s, f);
                sss.Width = Convert.ToInt32(size.Width);
                sss.Height = Convert.ToInt32(size.Height);
                gr.Dispose();
            }
            return sss;
        }

        public void AutoAdjustFont(double W)
        {
            _FontSize = W;
            float size = 1;
            do
            {
                defalutFont = new Font("Tahoma", size, FontStyle.Regular, GraphicsUnit.Point, ((byte)(0)));
                size++;
            } while (Math.Abs(getFontSize("●", defalutFont).Width - (W * Resolution)) > 5);
        }

        private void DrawBaseString(string s, float x, float y, Font f, Color color, bool Left = false)
        {
            Brush br = new SolidBrush(color);
            float dx = x * Resolution;
            float dy = y * Resolution;
            //計算文字長寬
            int img_width = getFontSize(s, f).Width, img_height = getFontSize(s, f).Height;
            if (Left)
            {
                g.DrawString(s, f, br, new PointF(GetRealX(dx), GetRealY(dy + img_height / 2)));
            }
            else
            {
                g.DrawString(s, f, br, new PointF(GetRealX(dx - img_width / 2), GetRealY(dy + img_height / 2)));

            }
        }

        private void DrawBaseLine(float x1, float y1, float x2, float y2, float w, Color color, bool IsDashLine = false)
        {
            float dx1 = x1 * Resolution;
            float dy1 = y1 * Resolution;
            float dx2 = x2 * Resolution;
            float dy2 = y2 * Resolution;
            float dw = w * Resolution;
            Brush br = new SolidBrush(color);
            Pen pen = new Pen(br, dw);
            if (IsDashLine)
                pen.DashStyle = DashStyle.DashDot;
            g.DrawLine(pen, GetRealX(dx1), GetRealY(dy1), GetRealX(dx2), GetRealY(dy2));


        }

        public void Dispose()
        {
            _Image.Dispose();
            _Image = null;
            GC.Collect();
        }
    }
}
