package StamFluidSolverPkg;

import processing.core.PApplet;

public class myBar extends myFluidObj {
	public float len;				//length of bar
	public myParticle A,B, Mp;		//endpoints of bar and midpoint
	public final float bLen = 150;			//default bar length
	public final float bMass = .003f;			//default bar mass
	public float bRad = 20f;			//default radius to scale against
	public float massRatio;				//percent of total mass in midpoint

	public myBar(StamFluidSolverGlob _p, float _x, float _y, float _mass, float _len) {
		super(_p, _x, _y, _mass);
		this.len = _len;
		A = new myParticle(_p, c.x-(.5f*len),c.y,bRad *.25f, _mass*.25f);
		B = new myParticle(_p, c.x+(.5f*len),c.y,bRad *.25f, _mass*.25f);
		Mp = new myParticle(_p, c.x, c.y, bRad *.5f, _mass*.5f);
		
		A.bounce = true;
		B.bounce = true;
		Mp.bounce = true;		
	}

	public myBar(StamFluidSolverGlob _p){
		this(_p, 0,0,0,0);
		len = bLen;
		mass = bMass;
		relocate(1);
		massRatio = .5f;
		distMassAndSize();
	}
		
	public void moveMe(float[] uAra, float[] vAra){
		//move each particle independently
		A.moveMe(uAra, vAra);
		B.moveMe(uAra, vAra);
		myVector mPold = new myVector(Mp.c);
		Mp.moveMe(uAra, vAra);
		if(p.drawCount != 1){
			myVector mPDiff = myVector._sub(Mp.c, mPold);
			A.U._add(Mp.U);
			B.U._add(Mp.U);
			A.c._add(mPDiff);
			B.c._add(mPDiff);
		}
		c.set(Mp.c);		
		//then resize line to maintain length
		reSizeLine();
	}

	public void drawMe(boolean setClr){
		//endpoints
		A.drawMe(setClr);
		B.drawMe(setClr);		
		Mp.drawMe(setClr);
		//line
		p.stroke(0,0,0);		//fade color?
		p.line(A.c.x,A.c.y,B.c.x,B.c.y);
		p.noStroke();
	}
	//override base class version so that mp is updated appropriately
	public void relocate(float mult){
		A.relocate(1.0f);
		B.relocate(1.0f);
		Mp.c = myVector._mult(myVector._add(A.c, B.c), .5f);
		reSizeLine();
	}
	
	public void reSizeLine(){//move A and B toward/away from midpoint to maintain length
		myVector AB = myVector._sub(B.c, A.c);//new myVector(B.c.x - A.c.x, B.c.y - A.c.y, 0);		//vec from A to B
		myVector newMp = myVector._add(myVector._mult(A.c, .5f), myVector._mult(B.c, .5f));
		Mp.setPos(newMp);
		c.set(Mp.c);			//make sure bar c follows Mp
		AB._normalize();
		AB._mult(.5f *len);
		myVector Ap = myVector._add(Mp.c, myVector._mult(AB, -1)),
				Bp = myVector._add(Mp.c, AB);
		A.setPos(Ap);
		B.setPos(Bp);	
		float dist = myVector._dist(Ap, Bp);
		if(PApplet.abs(len - dist) > p.epsVal){	
			System.out.println("len has changed : old "+ len + " new " + dist);			
		}
	}
	
	public void distMassAndSize(){//resdistribute mass to 3 component particles - 1/2 to midpoint, 1/4 to each end point
		float endRatio = (1.0f - massRatio) * .5f;
		A.mass = mass*endRatio;
		B.mass = mass*endRatio;
		Mp.mass = mass*massRatio;
		A.rad = bRad*endRatio;
		B.rad = bRad*endRatio;
		Mp.rad = bRad*massRatio;
	}
			
	public String toString(){		
		String result = super.toString(); 
		result += "bar length : " + len + " Endpoints A :" + A.toString() + " B : " + B.toString() + " Midpoint : " + Mp.toString();
		return result;	
	}
	
}
