package StamFluidSolverPkg;
/**
*  replacement class for pvector that is exportable for jni usage
*/
public class myVector{
  public float x,y,z;
  public float magn;
  
  myVector(float x, float y, float z){this.x = x; this.y = y; this.z = z; this.magn = (float) Math.sqrt((this.x*this.x) + (this.y*this.y) + (this.z*this.z));}         //constructor 3 args  
  myVector(myVector p){ this(p.x, p.y, p.z); }                                                                                                           //constructor 1 arg  
  myVector(){ this(0,0,0);}                                                                                                                               //constructor 0 args 
  public void set(float x, float y, float z){ this.x = x;  this.y = y;  this.z = z;  }                                                                     //set 3 args 
  public void set(myVector p){ this.x = p.x; this.y = p.y; this.z = p.z; }                                                                                //set 1 args
  public void _mult(float n){ this.x *= n; this.y *= n; this.z *= n;  }                                                                                    //_mult 3 args  
  public static myVector _mult(myVector p, float n){ myVector result = new myVector(p.x * n, p.y * n, p.z * n); return result;}                        //1 vec, 1 float
  public static myVector _mult(myVector p, myVector q){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); return result;}               //2 vec
  public static void _mult(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x *q.x, p.y * q.y, p.z * q.z); r.set(result);}       //2 vec src, 1 vec dest  
  public void _add(float x, float y, float z){ this.x += x; this.y += y; this.z += z;  }                                                                   //_add 3 args
  public void _add(myVector v){ this.x += v.x; this.y += v.y; this.z += v.z;  }                                                                           //_add 1 arg  
  public void _sub(float x, float y, float z){ this.x -= x; this.y -= y; this.z -= z;  }                                                                   //_sub 3 args
  public void _sub(myVector v){ this.x -= v.x; this.y -= v.y; this.z -= v.z;  }                                                                           //_sub 1 arg  
  public static myVector _add(myVector p, myVector q){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); return result;}                //2 vec
  public static myVector _sub(myVector p, myVector q){ myVector result = new myVector(p.x - q.x, p.y - q.y, p.z - q.z); return result;}                //2 vec
  public static void _add(myVector p, myVector q, myVector r){ myVector result = new myVector(p.x + q.x, p.y + q.y, p.z + q.z); r.set(result);}       //2 vec src, 1 vec dest  
  public float _mag(){ this.magn = (float) Math.sqrt((this.x*this.x) + (this.y*this.y) + (this.z*this.z)); return this.magn; }  
  public void _normalize(){this._mag();this.x /= this.magn; this.y /= this.magn; this.z /= this.magn;}
  
  public float _dist(myVector q){ return (float)Math.sqrt( ((this.x - q.x)*(this.x - q.x)) + ((this.y - q.y)*(this.y - q.y)) + ((this.z - q.z)*(this.z - q.z)) ); }
  public static float _dist(myVector q, myVector r){  return (float)Math.sqrt(((r.x - q.x) *(r.x - q.x)) + ((r.y - q.y) *(r.y - q.y)) + ((r.z - q.z) *(r.z - q.z)));}  
  public void _div(float q){this.x /= q; this.y /= q; this.z /= q;}  
  public myVector _cross(myVector b){ return new myVector((this.y * b.z) - (this.z*b.y), (this.z * b.x) - (this.x*b.z), (this.x * b.y) - (this.y*b.x));}//_cross product 
  public String toString(){return "(" + this.x + ", " + this.y + ", " + this.z+")";}
}//myPVect