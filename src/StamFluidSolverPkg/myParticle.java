package StamFluidSolverPkg;

/**
 * a particle to display in the fluid simulator
 * @author John
 */

public class myParticle extends myFluidObj {	
	public float rad;				//display radius of the particle
	public final float pMass = .001f;
	//radius of drawn particles
	public final float pRad = 3f;
	
	public boolean dying;
	public int lifespan;
	
	public myParticle(StamFluidSolverGlob _p, float _x, float _y, float _rad, float _mass) {
		super(_p,_x,_y, _mass);
		rad = _rad;		
	}	
	public myParticle(StamFluidSolverGlob _p, float _x, float _y, float _rad) {
		this(_p,_x,_y, _rad, 0);
		mass = pMass;		
	}	
	public myParticle(StamFluidSolverGlob _p){
		this(_p, 0,0,0,0);
		relocate(1);
		rad = pRad;
		mass = pMass;
	}
	
	public void setLifespan(int _lfspan){		dying = true;		lifespan = _lfspan;	}

	public void drawMe(boolean color){
		super.drawMe(color);
		p.pushMatrix();
			p.ellipse(c.x, c.y, rad, rad);
		p.popMatrix();	
	}
	
	public void setPos(myVector pos){c.x = pos.x; c.y = pos.y; c.z = pos.z;	}
	
	public String toString(){		
		String result = super.toString(); 
		result += "particle radius : " + rad;
		return result;	
	}
	
}//myParticle
