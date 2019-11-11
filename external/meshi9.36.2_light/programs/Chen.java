package programs;

import javax.swing.*;
import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: chen
 * Date: 09/01/14
 * Time: 21:37
 * To change this template use File | Settings | File Templates.
 */
public class Chen {

    public static void main(String[] argv) throws InstantiationException, IllegalAccessException {
        int[] toSort = {5, 4, 6 ,3 ,7, 2, 8, 1 ,9};
        int [] sorted = mergesort(toSort);
        for (int i:sorted)
            System.out.print(i+" ");
        System.out.println();

        Painter.draw("programs.Chen$Tree2");



    }
    static class Tree2 extends Tree {
        public Tree2() {
            super();
        }
        public Tree2(Pixel start, Pixel end) {
            super(start,end);
        }
       public void drawLine(){
           Painter.drawLine(start,end);
           if (start.getX()-end.getX()!= start.getY()-end.getY())
               Painter.drawLine(new Pixel(start.getX()+4,start.getY()+4),new Pixel(end.getX()+4,end.getY()+4));
           else
               Painter.drawLine(new Pixel(start.getX()+4,start.getY()),new Pixel(end.getX()+4,end.getY()));
       }
        public Tree getBranch(double theta) {
            Pixel newStart = end;
            Pixel tempPixel = end.translate(new Pixel(-start.getX(),-start.getY()));
            tempPixel = tempPixel.rotate(theta);
            tempPixel = new Pixel(tempPixel.getX()*0.5, tempPixel.getY()*0.5);
            Pixel newEnd = tempPixel.translate(end);
            return new Tree2(newStart,newEnd);
        }

    }
     static class Tree {
         protected Pixel start,end;
         double length;

         public Tree(Pixel start, Pixel end) {
             this.start = start;
             this.end = end;
             this.length = start.distanceFrom(end);
         }
         public Tree() {
             int width = Painter.getFrameHeight()/2;	// find the x coordinate of the center of the frame
             int height = Painter.getFrameWidth()/2; // find the y coordinate of the center of the frame

             start = new Pixel(width,0);
             end   = new Pixel(width,height);
             length = start.distanceFrom(end);
         }


        public void draw() {
            if (length < 1) return;
            Tree left = getBranch(-Math.PI/4);
            Tree right = getBranch(Math.PI/4);
            drawLine();
            left.draw();
            right.draw();
        }
         public Tree getBranch(double theta) {
             Pixel newStart = end;
             Pixel tempPixel = end.translate(new Pixel(-start.getX(),-start.getY()));
             tempPixel = tempPixel.rotate(theta);
             tempPixel = new Pixel(tempPixel.getX()*0.5, tempPixel.getY()*0.5);
             Pixel newEnd = tempPixel.translate(end);
             return new Tree(newStart,newEnd);
         }

         public void drawLine(){
             Painter.drawLine(start,end);
         }

    }
    public static int[] mergesort(int[] arr){
        if (arr.length<2)
            return arr;
        else{
            return
                    merge(mergesort(splitLeft(arr)),
                          mergesort(splitRight(arr)));
        }
    }
    public static void occurrences(int[] a){
            a = mergesort(a);
            int i = 0;
            while (i<a.length){
                int j = i+1;
                while (j<a.length && a[i]==a[j])
                    j = j+1;
                System.out.println(a[i] + " : " + (j-i));
                i=j;
            }
    }

    public static int[] splitLeft(int[] arr){
        int[] ans = new int[arr.length/2];
        for (int i=0; i<arr.length/2; i=i+1)
            ans[i]=arr[i];
        return ans;
    }
    public static int[] splitRight(int[] arr){
        int[] ans = new int[arr.length-arr.length/2];
        for (int i=arr.length/2; i<arr.length; i=i+1)
            ans[i-arr.length/2]=arr[i];
        return ans;
    }
    public static int[] merge(int[] arr1, int[] arr2) {
        int ind = 0, i1=0, i2=0;
        int len1 = arr1.length, len2=arr2.length;
        int[] ans = new int[len1+len2];
        while(i1<len1 & i2 < len2) {
            if(arr1[i1] <arr2[i2]) {
                ans[ind] = arr1[i1];
                i1=i1+1;
            }
            else {
                ans[ind] = arr2[i2];
                i2=i2+1;
            }
            ind=ind+1;
        }
        for(int i=i1;i<len1;i=i+1) {
            ans[ind] = arr1[i];
            ind=ind+1;
        }
        for(int i=i2;i<len2;i=i+1) {
            ans[ind] = arr2[i];
            ind=ind+1;
        }
        return ans;
    }


       private static class Painter extends JComponent {

        private static final long serialVersionUID = 1L;
        private static int SIZE = 600;
        private static Painter painter;
        private static Graphics graphics;
        private static String shape = null;


        // Create a frame and display it
        public static void draw(String aShape) {
            shape = aShape;
            JFrame frame = new JFrame(shape);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            frame.setLocationByPlatform(true);
            painter =  new Painter();
            frame.add(painter, BorderLayout.AFTER_LAST_LINE);
            frame.pack();
            frame.setVisible(true);
        }

        // returns the Frame's width
        public static int getFrameWidth () {
            return painter.getSize().width;
        }

        // returns the Frame's height
        public static int getFrameHeight () {
            return painter.getSize().height;
        }

        // changes the color of the lines to be drawn
        public static void setColor (String color) {
            if (color.equals("red")){
                graphics.setColor(Color.red);
            }
            else if (color.equals("blue")){
                graphics.setColor(Color.blue);
            }
            else if (color.equals("green")){
                graphics.setColor(Color.green);
            }
        }
           public static void drawLine (Pixel p1, Pixel p2) {
               drawLine((int)Math.round(p1.getX()),(int)Math.round(p1.getY()),(int)Math.round(p2.getX()),(int)Math.round(p2.getY()));

           }

           // Draw a line on the frame
           public static void drawLine (int x1, int y1, int x2, int y2) {
                graphics.drawLine(x1, getFrameHeight() - y1, x2, getFrameHeight() - y2);
           }

           // Set the default size of the window frame to SIZE*SIZE pixels
           public Dimension getPreferredSize() {
               return new Dimension(SIZE, SIZE);
           }

           // paint the frame - draw the shape given (call the draw method in that shape object)
           public void paint(Graphics g) {
               Painter.graphics = g;
               try{
                   Object myShape = (Class.forName(shape)).newInstance();
                   Object [] objs = null;
                   Class [] classes = null;
                   (Class.forName(shape)).getMethod("draw", classes).invoke(myShape, objs);
               }
               catch(Exception e)
               {
                   System.out.println("Can't handle shape " + shape);
                   System.out.println(e.toString());
                   System.out.println(e.getCause());

               }



           }
       }

        private static class Pixel {
            private double x;
            private double y;
            public static final int LEFT = -1;
            public static final int RIGHT = 1;
            public Pixel(){
                    x=0;
                    y=0;
            }
            public Pixel(double x, double y){
                this.x = x;
                this.y = y;
            }
            public double distanceFrom(Pixel other) {
                double dx = x - other.x;
                double dy = y - other.y;
                return Math.sqrt(dx*dx + dy*dy);
            }
            public double getX(){
                return x;
            }
            public double getY(){
                return y;
            }

            public Pixel rotateRelativeToPixel( Pixel p1,double theta){

                Pixel p3=this.translate( new Pixel(-p1.getX(),-p1.getY()));
                p3=p3.rotate(theta);
                p3=p3.translate( p1);
                return p3;
            }

            public Pixel translate(Pixel p){
                return new Pixel(x+p.x,y+p.y);
            }
            public Pixel rotate(double theta){
                return new Pixel(x*Math.cos(theta) - y*Math.sin(theta),
                                 x*Math.sin(theta) + y*Math.cos(theta));
            }
        }


    }

