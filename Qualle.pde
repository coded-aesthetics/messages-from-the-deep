import processing.opengl.*;
import ddf.minim.analysis.*;
import ddf.minim.*;


/**
 * What the wobble?! 
 * by Jan Reinsch based on the vertices example by Simon Greenwold.
 * What you see is a sphere with oscillating Gauss-functions
 * mapped on its surface, or is it?!
 */

class Vertex3D
{
  float x,y,z;
  float r, rho, phi;

  float lon, lat;

  Vertex3D(float x, float y, float z, boolean isCartesian)
  {
    if (isCartesian)
    {
      this.x = x;
      this.y = y;
      this.z = z;
      getSpheric();
    }
    else
    {
      this.r = x;
      this.rho = y;
      this.phi = z;

      getLonLat();      
      getCartesian();
    }
  }

  void setR(float r)
  {
    this.r = r;
    getCartesian();
  }

  void setRho(float rho)
  {
    this.rho = rho;
    getCartesian();
  }

  void setPhi(float phi)
  {
    this.phi = phi;
    getCartesian();
  }

  void getSpheric()
  {
    r = sqrt( x*x + y*y + z*z );
    if (y >= 0) rho = acos( x / sqrt( x*x + y*y ) );
    else rho = TWO_PI - acos( x / sqrt( x*x + y*y ) );
    phi = HALF_PI - atan( z / sqrt( x*x + y*y ) );

    getLonLat();
  }

  void getCartesian()
  {
    x = r * sin( phi ) * cos( rho );
    y = r * sin( phi ) * sin( rho );
    z = r * cos( phi );
  }

  void getLonLat()
  {
    this.lon = - ( TWO_PI - rho );
    this.lat = phi - HALF_PI;
  }

  float getSphericDist(Vertex3D v)
  {
    return r * acos(sin(v.lat)*sin(lat) + cos(v.lat)*cos(lat)*cos(lon - v.lon));
  }
}

class BeatListener implements AudioListener
{
  private BeatDetect beat;
  private AudioPlayer source;

  BeatListener(BeatDetect beat, AudioPlayer source)
  {
    this.source = source;
    this.source.addListener(this);
    this.beat = beat;
  }

  void samples(float[] samps)
  {
    beat.detect(source.mix);
  }

  void samples(float[] sampsL, float[] sampsR)
  {
    beat.detect(source.mix);
  }
}


class LFO
{
  float lon, amp, mid, phase;

  int timePassed;

  LFO(float lon, float amp, float mid, float phase)
  {
    this.lon = lon;
    this.amp = amp;
    this.mid = mid;
    this.phase = phase;

    timePassed = 0;
  }

  float tick(int ms)
  {
    float ret = sin(timePassed / lon * TWO_PI + phase) * amp + mid;
    timePassed += ms;
    return ret;
  }

  void reset(int timePassed)
  {
    this.timePassed = timePassed;
  }
}

class Color
{
  int red, green, blue, alpha;

  int col;

  Color(int col)
  {
    this.col = col;
    setColor(col);
  }

  Color(int r, int g, int b, int a)
  {
    red = r;
    green = g;
    blue = b;
    alpha = a; 

    this.col = 0x01000000 * alpha + 0x00010000 * red + 0x00000100 * green + blue;
  }

  void setColor(int col)
  {
    alpha  = (col & 0xFF000000) >> 24;
    red    = (col & 0x00FF0000) >> 16;
    green  = (col & 0x0000FF00) >>  8; 
    blue   = (col & 0x000000FF);
  }
}

class PerpetuumGradient
{
  Gradient gradOne, gradTwo;

  int lonOne, lonTwo, curLonOne, curLonTwo;

  PerpetuumGradient()
  {
    gradOne = new Gradient(round(random(0, 0xFFFFFF)) | 0xFF000000, 
    round(random(0, 0xFFFFFF)) | 0xFF000000);
    gradTwo = new Gradient(round(random(0, 0xFFFFFF)) | 0xFF000000, 
    round(random(0, 0xFFFFFF)) | 0xFF000000);   

    lonOne = round(random(1000, 10000));
    lonTwo = round(random(1000, 10000));  

    curLonOne = 0;
    curLonTwo = 0;
  }

  Gradient next(int ms)
  {
    if (curLonOne >= lonOne)
    {
      gradOne = new Gradient(gradOne.endColor.col, 
      round(random(0, 0xFFFFFF)) | 0xFF000000);
      curLonOne -= lonOne;
    }

    if (curLonTwo >= lonTwo)
    {
      gradTwo = new Gradient(gradTwo.endColor.col, 
      round(random(0, 0xFFFFFF)) | 0xFF000000);
      curLonTwo -= lonTwo;
    }

    curLonOne += ms;
    curLonTwo += ms;

    float amountOne = (float)curLonOne / (float) lonOne;
    float amountTwo = (float)curLonTwo / (float) lonTwo;

    return new Gradient(gradOne.getColor(amountOne), gradTwo.getColor(amountTwo));
  }
}

class Gradient
{
  Color startColor, endColor;

  Gradient (int startCol, int endCol)
  {
    setColor(startCol, endCol);
  }

  void setColor(int startCol, int endCol)
  {
    startColor = new Color(startCol);
    endColor   = new Color(endCol);
  }

  Gradient (Color startCol, Color endCol)
  {
    startColor = startCol;
    endColor   = endCol;
  }

  Color getColor(float amount)
  {
    int r = startColor.red   + round((float)(endColor.red   - startColor.red)   * amount);
    int g = startColor.green + round((float)(endColor.green - startColor.green) * amount);
    int b = startColor.blue  + round((float)(endColor.blue  - startColor.blue)  * amount);
    int a = startColor.alpha + round((float)(endColor.alpha - startColor.alpha) * amount);    

    return new Color(r, g, b, a);
  }
}

class EFunction
{
  float sigma, sc; 
  Vertex3D v;

  LFO lfoSC;

  EFunction(Vertex3D v, float sigma, float sc)
  {
    this.v = v;
    this.sigma = sigma;
    this.sc = sc;
  }

  void lfo(int lon, int amp, int mid, float phase)
  {
    this.lfoSC = new LFO(lon, amp, mid, phase);
  }

  void tick(int ms)
  {
    this.sc = lfoSC.tick(ms);
  }

  float getVal(Vertex3D v)
  {
    float distance = v.getSphericDist(this.v);

    if (distance > sigma*sigma) return 0;

    return sc * exp(-(distance*distance / (2.0*sigma*sigma)));
  }
}

/****************************************************************
 ****        VARIABLES // CONSTANTS      *************************
 ****************************************************************/

Minim minim;
AudioPlayer jingle;
FFT fft;
String windowName;

PImage img;

int NUM_E_FUNCTIONS = 15;

int NUM_FACES_HOR = 40;
int NUM_FACES_VER = 15;

float SCALE = 0.7;

float SPHERE_RADIUS = 200 * SCALE;

float MIN_E_SCALE = 50.0 * SCALE;
float MAX_E_SCALE = 100 * SCALE;

float MIN_E_SIGMA = 50.0 * SCALE;
float MAX_E_SIGMA = 75.0 * SCALE;

float count = 0;
float count_x = 0;
int lastTime = 0;
float particle_size = 1;

int cur_highlight_row = 0;


PImage particle_tex;

BeatDetect beat;
BeatListener bl;

PerpetuumGradient pg; // used to make colors change

Vertex3D [][] sphereVertices; // used to store the vertices of the sphere

EFunction [] e_funcs; // used to store all EFunctions used

LFO lfo, lfo2, lfo3;
Vertex3D lightPos;

boolean lightsOn = false;

float angY = 0;//5.8;
float angX = 0;//-2.44;

float offsetX = 0;
float offsetY = 0;
float lastMouseX = 0;
float lastMouseY = 0;

/****************************************************************
 ****        BEGIN OF FUNCTIONS          *************************
 ****************************************************************/

void mouseReleased()
{
  /*  lightsOn = !lightsOn;
   if (!lightsOn)
   {
   lastMouseX = mouseX - offsetX;//map(mouseY, 0, height, 0, -3*PI);
   lastMouseY = mouseY - offsetY;//map(mouseX, 0, width, 0, 3*PI);
   }
   else
   {
   offsetX = mouseX - lastMouseX;
   offsetY = mouseY - lastMouseY;
   }*/
}

void drawTexture()
{
  //int startColor = round(random(0, 0xFFFFFF)) | 0xFF000000;
  //int endColor   = round(random(0, 0xFFFFFF)) | 0xFF000000;  

  Gradient g = pg.next(50);
  LFO lfo = new LFO(120, 0.5, 0.5, 0);

  //img = loadImage("Foto208.jpg");
  img = createImage(120, 1, RGB);
  float amount = 0;
  for(int i=0; i < img.pixels.length; i++) {
    amount = lfo.tick(1);
    Color c = g.getColor(amount);//(float)i/(float)img.pixels.length);

    //Color c = g.getColor(amount);
    //print(c.red + " " + c.green + " " + c.blue + "//");
    img.pixels[i] = color(c.red, c.green, c.blue);//, 255);//i%img.width * 2); 
    //print (img.pixels[i] + "\n");
  }

  //img = loadImage("vintage_polka_dot_texture_by_aeiryn.jpg");
  //tint(255, 255, 255, 200);
}

void setup() {
  size(1290, 730, P3D);

  pg = new PerpetuumGradient();
  drawTexture();

  particle_tex = loadImage("particle.png");

  e_funcs = new EFunction[NUM_E_FUNCTIONS];

  for (int i = 0; i < NUM_E_FUNCTIONS; i++)
  {
    e_funcs[i] = new EFunction(new Vertex3D(SPHERE_RADIUS, random(0, TWO_PI), random(5*PI/12, 10*PI/12), false), random(MIN_E_SIGMA, MAX_E_SIGMA), random(35, 70)); 
    e_funcs[i].lfo(round(random(4000, 8000)), round(random(MIN_E_SCALE, MAX_E_SCALE)), 0, random(0, TWO_PI));
  }

  lastTime = millis();

  buildSphere(sphereVertices);

  lfo = new LFO(4200, 0.5, 0.5, 0);
  lfo2 = new LFO(3900*5, 0.5, PI, 0);
  lfo3 = new LFO(2300*5, 0.5, 0, 0);

  lightPos = new Vertex3D(300, 0, 0, false);

  minim = new Minim(this);

  jingle = minim.loadFile("soul implant.mp3", 2048);

  // create an FFT object that has a time-domain buffer the same size as jingle's sample buffer
  // note that this needs to be a power of two and that it means the size of the spectrum
  // will be 512. see the online tutorial for more info.
  fft = new FFT(jingle.bufferSize(), jingle.sampleRate());

  fft.linAverages(1);

  // a beat detection object that is FREQ_ENERGY mode that 
  // expects buffers the length of song's buffer size
  // and samples captured at songs's sample rate
  beat = new BeatDetect(jingle.bufferSize(), jingle.sampleRate());
  // set the sensitivity to 300 milliseconds
  // After a beat has been detected, the algorithm will wait for 300 milliseconds 
  // before allowing another beat to be reported. You can use this to dampen the 
  // algorithm if it is giving too many false-positives. The default value is 10, 
  // which is essentially no damping. If you try to set the sensitivity to a negative value, 
  // an error will be reported and it will be set to 10 instead. 
  beat.setSensitivity(300);
  //kickSize = snareSize = hatSize = 16;
  // make a new beat listener, so that we won't miss any buffers for the analysis
  bl = new BeatListener(beat, jingle);
  
  jingle.loop();
}

boolean go = true;

void draw() {
  //drawTexture();

  background(0);



  fft.forward(jingle.mix);

  //for(int i = 0; i < fft.specSize(); i++)
  //{
  // draw the line for frequency band i, scaling it by 4 so we can see it a bit better
  //line(i, height, i, height - fft.getBand(i)*4);
  //}
  //fill(255);
  //lights();
  //directionalLight(255, 255, 255, 0, -lfo.tick(millis() - lastTime), 0);

  //float fac = lfo.tick(millis() - lastTime);
  //float fac2 = lfo2.tick(millis() - lastTime);
  //float fac3 = lfo3.tick(millis() - lastTime);
  //float facSum = map(fac + fac2 + fac3, 0, 3.0, 0, 1);

  //directionalLight(100 + facSum * 155, 100 + facSum * 155, 100 + facSum * 155, 0, 0, -1);

  //if (lightsOn) directionalLight(255, 255, 255, 0, 0, -1);

  //pointLight(255, 255, 255
  //, -35, -fac * 2550 - 1280, -36);



  translate(width / 2, height / 2);

  angX = lfo2.tick(millis() - lastTime);
  angY = lfo3.tick(millis() - lastTime);  

  if (!lightsOn)
  {
    rotateY(angY);//map(mouseX, 0, width, 0, 3*PI));
    rotateX(angX);//map(mouseY, 0, height, 0, -3*PI));
  }
  else
  {
    angY = map(mouseX - offsetX, 0, width, 0, 3*PI);
    angX = map(mouseY - offsetY, 0, height, 0, -3*PI);
    rotateY(angY);
    rotateX(angX);
  }

  //println("angX: " + angX + " angY: " + angY + " mouseX: " + mouseX + " mouseY: " + mouseY);
  //println("rotate x: " + map(mouseY, 0, height, 0, -3*PI));  
  noStroke();
  fill(255, 255, 255);
  translate(0, 40, 100);

  if (true)
  {
    //lightPos.setRho(map(mouseX, 0, width, 0, 3*PI));
    //lightPos.setPhi(map(mouseY, 0, height, 0, -3*PI));
    lightPos.setRho(0);

    lightPos.setR(400);

    //pointLight(255, 255, 255
    //          , lightPos.x, lightPos.y, lightPos.z);
    lightPos.setRho(0);
    lightPos.setPhi(map(mouseY, 0, height, 0, -3*PI));  
    lightPos.setR(100);
    //pointLight(255, 255, 255
    // , lightPos.x, lightPos.y, lightPos.z);
    lightPos.setRho((count_x - count) * TWO_PI);
    lightPos.setPhi(count_x);  
    if ( beat.isKick() ) lightPos.setPhi(PI/4);
    lightPos.setR(100);
    pointLight(255, 255, 255
      , lightPos.x, lightPos.y, lightPos.z);

    //lightPos.setPhi(map(fft.getAvg(0), 0, 25, PI/4, 2*PI));
  }


  //if ( beat.isSnare() ) snareSize = 32;
  //if ( beat.isHat() ) hatSize = 32;

  SPHERE_RADIUS = map(fft.getAvg(0), 0, 10, 200 * SCALE, 280 * SCALE);
  //NUM_FACES_HOR = round(40 * map(fft.getAvg(0), 0, 25, 0, 1));

  for (int i = 0; i < NUM_E_FUNCTIONS; i++)
    e_funcs[i].tick(millis() - lastTime);

  lastTime = millis();


  addEFuncs();  

  drawSphere(SPHERE_RADIUS, NUM_FACES_HOR, NUM_FACES_VER);

  count += .007;

  count_x += .003;
}

void buildSphere(Vertex3D[][] vertices)
{
  float inc_hor = TWO_PI / (float)NUM_FACES_HOR;
  float inc_ver = PI / (float)NUM_FACES_VER;  
  float count_hor = 0;
  float count_ver = 0;

  Vertex3D curVertex, curVertex2;

  sphereVertices = new Vertex3D[NUM_FACES_VER + 2][];
  for (int i = 0; i < NUM_FACES_VER + 2; i++)
  {
    count_hor = 0;

    sphereVertices[i] = new Vertex3D[NUM_FACES_HOR + 1];
    for (int j = 0; j < NUM_FACES_HOR + 1; j++)
    {
      float rnd = 0;//random(-10, 10);

      sphereVertices[i][j] = new Vertex3D(SPHERE_RADIUS + rnd, count_hor, count_ver, false);
      count_hor += inc_hor;
    } 

    count_ver += inc_ver;
  }
}

void addEFuncs()
{

  buildSphere(sphereVertices);

  for (int i = 0; i < NUM_E_FUNCTIONS; i++)
    addE(e_funcs[i]);
}

void addE(Vertex3D v, float sigma, float sc)
{
  float inc_hor = TWO_PI / (float)NUM_FACES_HOR;
  float inc_ver = PI / (float)NUM_FACES_VER;  
  float ang_hor = 0;
  float ang_ver = 0;

  for (int i = round(NUM_FACES_VER / 2.5); i < NUM_FACES_VER + 1; i++)
  {
    ang_hor = 0; 

    for (int j = 0; j < NUM_FACES_HOR + 1; j++)
    {
      Vertex3D curVertex = new Vertex3D(SPHERE_RADIUS, ang_hor, ang_ver, false);

      float distance = curVertex.getSphericDist(v);
      println("dist: " + distance + " sigma: " + sigma);
      sphereVertices[i][j].setR(sphereVertices[i][j].r + sc * exp(-(distance*distance / (2.0*sigma*sigma))));
      //println (sphereVertices[i][j].r);
      ang_hor += inc_hor;
    }
    ang_ver += inc_ver;
  }
}

void addE(EFunction e)
{
  float inc_hor = TWO_PI / (float)NUM_FACES_HOR;
  float inc_ver = PI / (float)NUM_FACES_VER;  
  float ang_hor = 0;
  float ang_ver = 0;

  for (int i = 0; i < NUM_FACES_VER + 1; i++)
  {
    ang_hor = 0; 

    for (int j = 0; j < NUM_FACES_HOR + 1; j++)
    {
      Vertex3D curVertex = new Vertex3D(SPHERE_RADIUS, ang_hor, ang_ver, false);

      //float h = e.getVal(curVertex);

      sphereVertices[i][j].setR(sphereVertices[i][j].r + e.getVal(curVertex));
      //println (sphereVertices[i][j].r);
      ang_hor += inc_hor;
    }
    ang_ver += inc_ver;
  }
}

void drawSphere(float r, int faces_hor, int faces_ver)
{
  float inc_hor = TWO_PI / (float)faces_hor;
  float inc_ver = PI / (float)faces_ver;  
  float count_hor = 0;
  float count_ver = 0;

  Vertex3D curVertex, curVertex2;

  if (beat.isSnare()) particle_size = random(0.5, 1);

  for (int i = round(faces_ver / 2.5); i < faces_ver; i++)
  {    
    //beginShape(QUAD_STRIP);
    //texture(img);
    //count_hor = 0;

    for (int s = 0; s < faces_hor + 1; s++) 
    {
      curVertex  = sphereVertices[i]    [s];
      curVertex2 = sphereVertices[i + 1][s];

      if (curVertex.r < 0) curVertex.setR(0);
      if (curVertex2.r < 0) curVertex2.setR(0);
      pushMatrix();
      translate(-curVertex.x, -curVertex.y, -curVertex.z);
      //rotateX(random(0, TWO_PI));

      //popMatrix();
      rotateX(-angX);
      rotateY(-angY);

      //float fac = map(mouseX, 0, width, faces_ver - i - 1, 1);
      float fac = map(particle_size, 0, 1, faces_ver - i - 1, 1);// random(0.5, 1);
      //println(fac);
      //popMatrix();
      beginShape(TRIANGLE_STRIP);
      texture(particle_tex);

      float bla = map(i, round(faces_ver / 2.5), faces_ver, 0, 255);
      if ((s + cur_highlight_row) % 5 == 0) bla = 255;
      tint(255, bla);
      vertex(-15 / fac, -15/fac, 0,  0,  0);
      vertex(+15/fac, -15/fac, 0, 30,  0);
      vertex(-15/fac, +15/fac, 0,  0, 30);

      vertex(+15/fac, +15/fac, 0, 30, 30);

      //translate(curVertex.x, curVertex.y, curVertex.z);
      endShape();   
      translate(curVertex.x, curVertex.y, curVertex.z);
      popMatrix();

      //float texCoordX = curVertex.r / (float)(SPHERE_RADIUS * 2.6);
      //if (texCoordX >= 0.97) texCoordX = 0.97;

      //float texCoordX2 = curVertex2.r / (float)(SPHERE_RADIUS * 2.6);
      //if (texCoordX2 >= 0.97) texCoordX2 = 0.97;

      //textureMode(NORMALIZED);
      //vertex(curVertex.x,  curVertex.y,  curVertex.z, (float) i / (float) (faces_ver +1) * img.width, (float)s / (float) (faces_hor + 1) * img.height);
      //vertex(curVertex2.x, curVertex2.y, curVertex2.z, (float) (i+1) / (float) (faces_ver +1) * img.width, (float)s / (float) (faces_hor + 1) * img.height);       
      //vertex(curVertex.x,  curVertex.y,  curVertex.z, texCoordX, 0.5);
      //vertex(curVertex2.x, curVertex2.y, curVertex2.z, texCoordX2, 0.5);       


      count_hor += inc_hor;
    }

    curVertex  = sphereVertices[i]    [cur_highlight_row];
    curVertex2 = sphereVertices[i + 1][cur_highlight_row];

    vertex(curVertex.x,  curVertex.y,  curVertex.z);//, (float) i / (float) (faces_ver + 1) * img.width, img.height);
    vertex(curVertex2.x, curVertex2.y, curVertex2.z);//, (float) (i+1) / (float) (faces_ver + 1) * img.width, img.height);       

    endShape();

    count_ver += inc_ver;
  }
  if (beat.isHat()) cur_highlight_row = (cur_highlight_row + 1) % 5;
}

