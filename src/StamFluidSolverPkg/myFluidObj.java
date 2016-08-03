package StamFluidSolverPkg;

import processing.core.PApplet;
import processing.core.PConstants;

public abstract class myFluidObj {
	public StamFluidSolverGlob p;
	public int ID;	
	public static int objCount = 0;
	public float mass;
	public int r,g,b;				
	public myVector c;				//object's COM loc
	public myVector U;				//object's COM vel
	
	public boolean relObj;			//relocate obj COM to some random central screen position if boundary is hit
	public boolean bounce;			//if should bounce away from boundary
	
	//multipliers to change the color of an object as it moves
	public final float objVelClrMultX = 5000;
	public final float objVelClrMultY = 5000;	
	
	
	public myFluidObj(StamFluidSolverGlob _p, float _x, float _y, float _mass) {
		ID = objCount++;
		p = _p;
		c = new myVector(_x,_y,0);
		mass = _mass;
		U = new myVector(0,0,0);	
		r = 0; g = 0; b = 0;
	}	
	
	public float calcRandLocation(float randNum1, float randNum2, float sketch, float mathCalc, float mult){
		float retVal = (sketch/2.0f) + (randNum2 * (sketch/3.0f) * mathCalc * mult);
		return retVal;
	}

	/**
	 * this will subject this object to fluid forces
	 */
	public void moveMe(float[] uAra, float[] vAra){
		float tmpX = c.x, tmpY = c.y;
		
		if((tmpX < 0) || (tmpX >= p.wLessCSize)){			tmpX = (int)(tmpX < 0 ? 0 : p.wLessCSize);		}
		if((tmpY < 0) || (tmpY >= p.hLessCSize)){			tmpY = (int)(tmpY < 0 ? 0 : p.hLessCSize);		}				
		
		int cellXIDX = (int)((tmpX)/p.cellWidth), cellYIDX = (int)((tmpY)/p.cellHeight);
		U.x = uAra[p.IX(cellXIDX, cellYIDX)];		U.y = vAra[p.IX(cellXIDX, cellYIDX)];		
		tmpX += U.x * p.deltaT / mass;				tmpY += U.y * p.deltaT / mass;
		
		if((tmpX <= 0) || (tmpX >= p.wLessCSize)){	
			if(relObj){relocate(.5f);} 
			else {
				if(tmpX <= 0){		c.x = -tmpX;	U.x = (bounce? U.x * (tmpX - 1.0f) : 0);	
				} else {
					float diff = p.wLessCSize - tmpX;//less than 0
					c.x = p.wLessCSize + diff;		U.x = (bounce? U.x * (diff - 1.0f) : 0);	
				}		
			}
		} else {			c.x = tmpX;		}
		if((tmpY <= 0) || (tmpY >= p.hLessCSize)){	
			if(relObj){relocate(.5f);}
			else {
				if(tmpY <= 0) {		c.y = -tmpY;	U.y = (bounce? U.y * (tmpY - 1.0f) : 0);		//reverses direction
				} else {
					float diff = p.hLessCSize - tmpY;//less than 0
					c.y = p.hLessCSize + diff;		U.y = (bounce? U.y * (diff - 1.0f) : 0);	
				}
			}
		} else {			c.y = tmpY;		}
	}//moveMe		
	//move object to random location in circle around center of screen
	public void relocate(float mult){
		float randNum1 = p.random(1);
		float randNum2 = p.random(1);
		c.x = calcRandLocation(randNum1, randNum2, p.sketchX, PApplet.sin(2 * PConstants.PI * randNum1), mult );
		c.y = calcRandLocation(randNum1, randNum2, p.sketchY, PApplet.cos(2 * PConstants.PI * randNum1), mult );	
	}
	//changes color of object depending on velocity
	public void colorMe(){
		r = (int) (PApplet.abs(U.x) * objVelClrMultX); 
		g = (int) (PApplet.abs(U.y) * objVelClrMultY);			
	}
	
	public void drawMe(boolean color){
		if(color){colorMe();}
		p.fill(r,g,b);
	}
	/**
	 * this will handle verifying an object's behavior at a boundary defined by 0 and w - bind it to boundary coord if it exceeds boundary
	 * @param x coord
	 * @param w large boundary
	 * @return modified coord value based on boundary
	 */
	public float handleBounds(float x, int w){return ((x <= 0) || (x >= w) ? -1 : x );}	

	public String toString(){
		String result = "Object : " + ID + " coords (" + c.x + "," + c.y +")";
		return result;	
	}

}
