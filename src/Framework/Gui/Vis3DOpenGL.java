//package Framework.Gui;
//
//import org.lwjgl.LWJGLException;
//import org.lwjgl.opengl.Display;
//import org.lwjgl.opengl.DisplayMode;
//
//import static org.lwjgl.opengl.GL11.*;
//
///**
// * Created by rafael on 5/28/17.
// */
//public class Vis3DOpenGL{
//    final boolean active;
//    public final int xDim;
//    public final int yDim;
//    public final int zDim;
//    final int maxDim;
//    final float scaleDims;
//    final float transXY;
//    final float transZ;
//    public final int xPix;
//    public final int yPix;
//
//    public Vis3DOpenGL(int xPix, int yPix, int xDim, int yDim, int zDim, String title, boolean active) {
//        this.active = active;
//        int maxDim = Math.max(xDim, yDim);
//        this.maxDim = Math.max(maxDim, zDim);
//        this.xDim = xDim;
//        this.yDim = yDim;
//        this.zDim = zDim;
//        this.xPix = xPix;
//        this.yPix = yPix;
//        scaleDims = (float) (2.0 / this.maxDim);
//        transXY = (float) (-this.maxDim / 2.0);
//        transZ = (float) (-zDim * 0.6);
//
//        if (active) {
//            try {
//                Display.setDisplayMode(new DisplayMode(xPix, yPix));
//                Display.setTitle(title);
//                Display.create();
//            } catch (LWJGLException e) {
//                e.printStackTrace();
//                System.err.println("unable to create Vis3D display");
//            }
//            glEnable(GL_DEPTH_TEST);
//            glMatrixMode(GL_PROJECTION);
//            glLoadIdentity();
//            glFrustum(-1, 1, -1, 1, 1, 3);
//            //glFrustum(0,maxDim,0,maxDim,maxDim,maxDim+zDim);
//            glMatrixMode(GL_MODELVIEW);
//            glLoadIdentity();
//            glScalef(scaleDims, scaleDims, scaleDims);
//            glTranslatef(transXY, transXY, transZ);
//        }
//    }
//    public Vis3DOpenGL(int xPix, int yPix, int xDim, int yDim, int zDim, String title) {
//        this.active = true;
//        int maxDim = Math.max(xDim, yDim);
//        this.maxDim = Math.max(maxDim, zDim);
//        this.xDim = xDim;
//        this.yDim = yDim;
//        this.zDim = zDim;
//        this.xPix = xPix;
//        this.yPix = yPix;
//        scaleDims = (float) (2.0 / this.maxDim);
//        transXY = (float) (-this.maxDim / 2.0);
//        transZ = (float) (-zDim * 0.6);
//
//        if (active) {
//            try {
//                Display.setDisplayMode(new DisplayMode(xPix, yPix));
//                Display.setTitle(title);
//                Display.create();
//            } catch (LWJGLException e) {
//                e.printStackTrace();
//                System.err.println("unable to create Vis3D display");
//            }
//            glEnable(GL_DEPTH_TEST);
//            glMatrixMode(GL_PROJECTION);
//            glLoadIdentity();
//            glFrustum(-1, 1, -1, 1, 1, 3);
//            //glFrustum(0,maxDim,0,maxDim,maxDim,maxDim+zDim);
//            glMatrixMode(GL_MODELVIEW);
//            glLoadIdentity();
//            glScalef(scaleDims, scaleDims, scaleDims);
//            glTranslatef(transXY, transXY, transZ);
//        }
//    }
//    public void Clear(float r,float g,float b){
//        if(active) {
//            glClearColor(r, g, b, 1);
//            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//        }
//    }
//    public void Show(){
//        if(active) {
//            Display.update();
//        }
//    }
//    public boolean CheckClosed(){
//        if(active) {
//            return Display.isCloseRequested();
//        }
//        return true;
//    }
//    public void Dispose(){
//        if(active) {
//            Display.destroy();
//        }
//    }
//    public boolean IsActive(){
//        return active;
//    }
//    public void FanShape(float centerX,float centerY,float centerZ,float scale,float[]points,float r,float g,float b) {
//        if(active) {
//            glColor3f(r, g, b);
//
//            //glLoadIdentity();
//            //glScalef(scaleDims,scaleDims,scaleDims);
//            //glScalef((float)(2*maxDim),(float)(2*maxDim),(float)(zDim/2.0));
//            //glScalef((float)(2*maxDim),(float)(2*maxDim),(float)(zDim/2.0));
//            //glTranslatef(-2,-2,-5);
//            //glTranslatef(transXY,transXY,transZ);
//            glBegin(GL_TRIANGLE_FAN);
//            glVertex3f((centerX), (centerY), -(centerZ));
//            for (int i = 0; i < points.length / 2; i++) {
//                float x = ((points[i * 2] * scale + centerX));
//                float y = ((points[i * 2 + 1] * scale + centerY));
//                float z = -((centerZ));
//                glVertex3f(x, y, z);
//            }
//            glEnd();
////        glColor3f(r,g,b);
////        glBegin(GL_TRIANGLE_FAN);
////        glLoadIdentity();
////        glTranslatef(-1,-1,0);
////        glRotatef(rotation,0,1,0);
////        glTranslatef(0,0,2);
////        glVertex3f(centerX,centerY,centerZ);
////        for (int i = 0; i < points.length/2 ; i++) {
////            float x=points[i*2]*scale+centerX;
////            float y=points[i*2+1]*scale+centerY;
////            float z=centerZ;
////            glVertex3f(x,y,z);
////        }
////        glEnd();
//        }
//    }
//
//    public void FanShapeOld(float centerX,float centerY,float centerZ,float scale,float[]points,float r,float g,float b) {
//        if (active) {
//            glColor3f(r, g, b);
//            glBegin(GL_TRIANGLE_FAN);
//            glVertex3f((centerX / xDim) * 2 - 1, (centerY / yDim) * 2 - 1, -((centerZ / zDim) + 2));
//            for (int i = 0; i < points.length / 2; i++) {
//                float x = ((points[i * 2] * scale + centerX) / xDim) * 2 - 1;
//                float y = ((points[i * 2 + 1] * scale + centerY) / yDim) * 2 - 1;
//                float z = -((centerZ / zDim) + 2);
//                glVertex3f(x, y, z);
//            }
//            glEnd();
//        }
//    }
//}
