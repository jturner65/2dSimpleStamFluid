package StamFluidSolverPkg;

import processing.core.*;

public class StamFluidSolverGlob extends PApplet {
	/**
	 * global variables and functions
	 * 
	 */
	// boolean flags used to control various elements of the program
	private boolean[] flags;
	// whether the sim should run or not
	public final int runSim = 0;
	// whether to display particles
	public final int dispParticles = 1;
	// whether or not to save an animation
	public final int saveAnim = 2;
	// performing concentric counter-force revolutions
	public final int concForce = 3;
	// eject particles and force from mouse pointer like magic wand
	public final int magicWand = 4;
	// number of flags in the boolean flags array
	private int numFlags = 5;
	// sim names for flags values output function
	public final String[] flagNames = { "runSim", "dispParticles", "saveAnim","saturn sim", "magic wand" };

	// counter for draw cycles
	public int drawCount;
	// how many cycles before a draw
	private int cycleModDraw;
	//radius for saturn sim calc
	float radsq;
	
	// path and filename to save pictures for animation
	public String animPath;
	public String animFileName;
	public int animCounter;

	// delta t increment value
	public float deltaT = .1f;

	//stam fluid solver
	public final float epsVal = .0001f;
	//arrays for density, xVelocity, yVelocity, now and for previous time step(for advection)
	public float[]  densityAra, oldDensityAra, uAra, uOldAra, vAra, vOldAra;
	public int[] cellCentersX, cellCentersY;
	public myFluidObj[] objAra;

	// size of cells and sketch
	public int cellSize;
	public int sketchX;
	public int sketchY;
	
	public int wLessCSize, hLessCSize;
	//size of array
	public int araSize;
	//number of cells across x and along y, leaving a border around the edges for boundaries
	public int numCellsX;
	public int numCellsY;
	//diffusion rate for fluid
	public final float diffRate = 0.0012f;
	//viscosity value for stam solver
	public float visc = .001f;
		
	//introduced force by clicking on screen
	public final float msDragForce = 7f;	
	//max force from clicking magicwand style
	public final float clickMaxForce = 500f;
	//force incrementer
	public final float fcIncr = 5f;
	//force from click
	public float clickForce;
	
	//number of particles in scene
	public final int numParticles = 5000;
	//number of lines in scene
	public final int numBars = 3;	
	//density "strength" at click
	public final float densityStr = 50;
	//radius of density
	public final int densityRadius = 4;
	//multiplier for density to fill color value
	public float densityColorMult;
	//for density calculation
	public float cellWidth;
	public float cellHeight;

	public int msCalcX;
	public int msCalcY;
	
	public void settings(){	size(800, 800, P3D); }// with cells of 4x4 pxls, this is 200 x 200 cells}
	public void setup() {
	//	size(800, 800, OPENGL); // with cells of 4x4 pxls, this is 200 x 200 cells
		this.frameRate(60.0f);
		this.colorMode(RGB, 255);
		this.background(255);
		this.noStroke();
		this.sketchX = 800;
		this.sketchY = 800;
		initOnce();
	}// setup

	/**
	 * draw
	 */
	// Draw the scene
	public void draw() {
		background(255);
		// initDisplay();
		if (flags[runSim]) {
			calcStamCycle();
		}
		drawCount++;
		if (drawCount % cycleModDraw == 0) {
			drawStamCycle();
			if (flags[saveAnim]) {
				save(animPath + animFileName + ((animCounter < 10) ? "000" : ((animCounter < 100) ? "00" : ((animCounter < 1000) ? "0" : ""))) + animCounter + ".jpg");
				animCounter++;
			}
			//drawCount = 0;
		}// if drawcount mod cyclemoddraw = 0
	}// draw

	/**
	 * initialize the program on setup
	 */
	private void initOnce() {
		flags = new boolean[numFlags];
		for (int i = 0; i < this.numFlags; ++i) { flags[i] = false; }

		flags[dispParticles] = true;
		// set sim to run initially
		flags[runSim] = true;
		flags[magicWand] = true;
		cellSize = 4;
		this.numCellsX = ((this.sketchX)/this.cellSize) - 2;
		radsq = (this.numCellsX/4.0f) * (this.numCellsX/4.0f);				
		this.numCellsY = ((this.sketchY)/this.cellSize) - 2;	
		this.araSize = (this.numCellsX + 2) * (this.numCellsY + 2);
		
		wLessCSize = (width - cellSize);
		hLessCSize = (height - cellSize);
		
		initProgram();
	}// initOnce
	
	/**
	 * reinitialize program due to user input
	 */
	public void initProgram() {
		drawCount = 0;
		cycleModDraw = 1;
		
		//stam solver fields
			//center of each cell, to draw flow lines
		this.cellCentersX 	= new int[this.araSize];
		this.cellCentersY 	= new int[this.araSize];
			//determine the centers of the cells for each cell in grid
		for(int i = 0; i < this.numCellsX+2; ++i){
			for(int j = 0; j < this.numCellsY+2; ++j){
				this.cellCentersX[IX(i,j)] = (int)( (i + 0.5f) * this.cellSize);
				this.cellCentersY[IX(i,j)] = (int)( (j + 0.5f) * this.cellSize);
			}
		}
		
		this.densityAra 	= new float[this.araSize];
		this.oldDensityAra 	= new float[this.araSize];
		this.uAra 			= new float[this.araSize];
		this.uOldAra 		= new float[this.araSize];
		this.vAra 			= new float[this.araSize];
		this.vOldAra 		= new float[this.araSize];
		this.cellWidth = (1.0f*this.width) / (this.numCellsX + 2.0f);
		this.cellHeight = (1.0f*this.height) / (this.numCellsY + 2.0f);
		
		this.densityColorMult = 40f;

		int numObjs = this.numParticles + this.numBars;
		this.objAra = new myFluidObj[numObjs];
		for(int i = 0; i < this.numParticles; ++i){ this.objAra[i] = new myParticle(this);}
		for(int i = this.numParticles; i < numObjs; ++i){ this.objAra[i] = new myBar(this);}
						
		animPath = sketchPath() + "\\" + (flags[dispParticles] ? "partAnim" + (int) random(10) : "flowAnim" + (int) random(10));
		animFileName = "\\" + "img" + (flags[dispParticles] ? "partAnim"	: "flowAnim");
		animCounter = 0;

		showFlags();
	}// initProgram
	
	/**
	 * calls a single cycle of the stam fluid solver
	 */
	public void calcStamCycle(){
		this.queryUI(this.oldDensityAra, this.uOldAra, this.vOldAra);
		this.calcVelocityStep(this.uAra, this.vAra, this.uOldAra, this.vOldAra);
		this.calcDensityStep(this.densityAra, this.oldDensityAra, this.uAra, this.vAra);
		this.moveObjects();
	}
	
	/**
	 * draws the results of the stam fluid solver
	 */
	public void drawStamCycle(){
		//this.drawDensity();
		if (flags[dispParticles]){		this.drawFlowLines(false);  this.drawObjects();} 
		else {							this.drawFlowLines(true);}
	}
	
	/**
	 * query ui for any mouse input
	 * @param oldDensityAra densities at last time step
	 * @param uOldAra x-dir velocity at last time step
	 * @param vOldAra y-dir velocity at last time step
	 */
	public void queryUI(float[] oldDensityAra, float[] uOldAra, float[] vOldAra){
		for (int i = 0; i < this.araSize; ++i){
			uOldAra[i] = 0f;
			vOldAra[i] = 0f;
			oldDensityAra[i] = 0f;
		}
		if (mousePressed){
			//find cell clicked in
			int calcX = (int)( min(1.0f, (mouseX / ((float)width-1)))  * this.numCellsX );
			int calcY = (int)( min(1.0f, (mouseY / ((float)height-1))) * this.numCellsY );
			
//			int calcX = (int)( (mouseX / (float)this.width)  * this.numCellsX );
//			int calcY = (int)( (mouseY / (float)this.height) * this.numCellsY );
			calcX = ( calcX < 1 || calcX > this.numCellsX + 1) ? ((calcX < 1 ) ? 1 : this.numCellsX+1) : calcX;
			calcY = ( calcY < 1 || calcY > this.numCellsY + 1) ? ((calcY < 1 ) ? 1 : this.numCellsY+1) : calcY;
			float sgn = 0;
				
			//if ( mouseButton == RIGHT ) { this.calcDensityEffect(calcX, calcY, oldDensityAra);}	
			if ( mouseButton == LEFT ) {		sgn = 1.0f;}
			if ( mouseButton == RIGHT ) { 		sgn = -1.0f;}
			uOldAra[IX(calcX,calcY)] = this.msDragForce * sgn * (mouseX - pmouseX);
			vOldAra[IX(calcX,calcY)] = this.msDragForce * sgn * (mouseY - pmouseY);
			
			//if(flags[magicWand]){
				//clickForce = (clickForce >= this.clickMaxForce ? this.clickMaxForce : clickForce + fcIncr);
			    //float val = sgn * clickForce;
			    //apply force in x and y to particles directly, scaled by dist
			//}//
			
		}//if mousePressed
		else {
			clickForce = 0;			//no click force when no click
			if (flags[concForce]){addCounterForces(uOldAra, vOldAra); }			//calc concentric forces
		}
	}//process gui input during sim
	
	//attempt to model hexagon on saturn by adding forces out to specific radius (in cells) in cw dir, and past radius in ccw dir
	public void addCounterForces(float[] uOldAra, float[] vOldAra){
		for (int x = 0; x < this.numCellsX; ++x){
			for (int y = 0; y < this.numCellsY; ++y){
				float thisSqRad = (((x - this.numCellsX/2) * (x - this.numCellsX/2)) +  ((y - this.numCellsY/2) * (y - this.numCellsY/2)));
				//float scaleFact = thisSqRad/radsq;
				if((radsq *.8f < thisSqRad) && (radsq * 1.2f >  thisSqRad)) {//ring circle
					uOldAra[IX(x,y)] = -(y - this.numCellsY/2) * 18f/radsq;//x direction use y
					vOldAra[IX(x,y)] = (x - this.numCellsX/2) * 18f/radsq;//y direction use x
				}//if in circle
				else if(radsq *.8f > thisSqRad) {//inner circles
					uOldAra[IX(x,y)] = ((y - this.numCellsY/2)) * 3f/radsq;//x direction use y
					vOldAra[IX(x,y)] = (-(x - this.numCellsX/2)) * 3f/radsq;//y direction use x
				}
				else{//big circles
					uOldAra[IX(x,y)] = ((y - this.numCellsY/2) - .2f*(x - this.numCellsX/2)) * 3f/radsq;//x direction use y
					vOldAra[IX(x,y)] = (-(x - this.numCellsX/2) - .2f*(y - this.numCellsY/2)) * 3f/radsq;//y direction use x
				}
			}//for y			
		}//for x
	}//addCounterForces
	
	/**
	 * this will move objects around in field
	 */
	public void moveObjects(){
		for(myFluidObj tmpObj : this.objAra){ 
			tmpObj.moveMe(this.uAra, this.vAra);
		}	
	}
	
	
	/**
	 * this will draw the particles and bars moving around in the fluid
	 */
	public void drawObjects(){
		pushMatrix();
			for(myFluidObj tmpObj : this.objAra){ 
				//tmpObj.moveMe(this.uAra, this.vAra);
				tmpObj.drawMe(!flags[concForce]);
			}
		popMatrix();
	}//drawObjects
	
	/**
	 * this will draw a physical representation of the velocity vector at each cell
	 */
	public void drawFlowLines(boolean dark){
		pushMatrix();
			stroke(dark ? 255 : 220);
			strokeWeight(1);
			float x0, y0, x1, y1;
			fill(dark ? 0 : 255);
			for(int i = 0; i < this.numCellsX+1; ++i){
				for(int j = 0; j < this.numCellsY+1; ++j){
					x0 = this.cellCentersX[IX(i,j)];
					y0 = this.cellCentersY[IX(i,j)];
					x1 = x0 + (4*uAra[IX(i,j)]*(cellSize+2)/deltaT);//((abs(4*uAra[IX(i,j)]*(cellSize+2)/deltaT) > epsVal) ? (4*uAra[IX(i,j)]*(cellSize+2)/deltaT) : 0 );
					y1 = y0 + (4*vAra[IX(i,j)]*(cellSize+2)/deltaT);//((abs(4*vAra[IX(i,j)]*(cellSize+2)/deltaT) > epsVal) ? (4*vAra[IX(i,j)]*(cellSize+2)/deltaT) : 0 );
					line(x0, y0, x1, y1);
				}//for j
			}//for i	
			noStroke();
		popMatrix();
	}//drawFlowLines
	
	/**
	 * will draw the fluid density result from clicking
	 */
	public void drawDensity(){
		for(int i = 0; i < this.numCellsX+2; ++i){
		    for(int j = 0; j < this.numCellsY+2; ++j){
		    	float densityVal = this.densityAra[IX(i,j)] * densityColorMult; 
		    	if(flags[dispParticles]){
		    		fill(255-((densityVal >255 ) ? (255) : (densityVal)  ), 255- ((densityVal >255 ) ? (255) : (densityVal)  ), 255);
		    	} else {
		    		fill(0,0,((densityVal >255 ) ? (255) : (densityVal)  ));		    		
		    	}
		    	rect(this.cellWidth * i, this.cellHeight * j , this.cellWidth, this.cellHeight);
		    }//for j
		}//for i
	}//drawDensity
	
	//jos stam solver functions	
	/**
	 * calculates the diffusion rate using a implicit integeration with gauss seidel relaxation
	 * @param bound boundary variable
	 * @param xAra array of densities/velocities
	 * @param xOldAra array of old densities/velocities
	 * @param diff diffusion rate
	 * @param dt
	 */
	public void calcDiffusionJS(int bound, float[] xAra, float[] xOldAra, float diffRate ){
		float a = this.deltaT * diffRate * this.numCellsX * this.numCellsY;
		
		for (int k=0 ; k<20 ; k++ ) {
			for (int i=1 ; i <= this.numCellsX; i++ ) {
				for (int j = 1 ; j <= this.numCellsY; j++ ) {
					xAra[IX(i,j)] = (xOldAra[IX(i,j)] + a * (xAra[IX(i-1,j)] + xAra[IX(i+1,j)] + xAra[IX(i,j-1)] + xAra[IX(i,j+1)]))/(1 + (4.0f * a));
				}//for j
			}//for i
			setBoundariesJS(bound, xAra );
		}//for k
	}//calcDiffuse
	
	/**
	 * handles boundary conditions for values held in an array
	 * @param bound particular boundary we are handling 
	 * @param xAra array of values being processed
	 */
	public void setBoundariesJS( int bound, float[] xAra ){	
		for (int  idx=1 ; idx < this.numCellsX  ; ++idx ) {
			xAra[IX(0, idx)] 				= ((bound==1) ? (-xAra[IX(1,idx)]) 				: (xAra[IX(1,idx)]));
			xAra[IX(this.numCellsX+1, idx)] = ((bound==1) ? (-xAra[IX(this.numCellsX,idx)]) : (xAra[IX(this.numCellsX,idx)]));
			xAra[IX(idx, 0 )] 				= ((bound==2) ? (-xAra[IX(idx,1)]) 				: (xAra[IX(idx,1)]));
			xAra[IX(idx, this.numCellsX+1)] = ((bound==2) ? (-xAra[IX(idx,this.numCellsX)]) : (xAra[IX(idx,this.numCellsX)]));
		}
		xAra[IX(0, 0 )]									 		= 0.5f * (xAra[IX(1,0 )] + xAra[IX(0 ,1)]);
		xAra[IX(0, this.numCellsX + 1)] 						= 0.5f * (xAra[IX(1,this.numCellsX + 1)] +xAra[IX(0 ,this.numCellsX )]);
		xAra[IX(this.numCellsX + 1, 0 )] 						= 0.5f * (xAra[IX(this.numCellsX, 0 )]+xAra[IX(this.numCellsX + 1, 1)]);
		xAra[IX(this.numCellsX + 1, this.numCellsX + 1)] 		= 0.5f * (xAra[IX(this.numCellsX,this.numCellsX + 1)] + xAra[IX(this.numCellsX + 1, this.numCellsX )]);
	}//setBoundaries
	
	/**
	 * calculates the advection step first order method
	 * @param bound boundary conditions
	 * @param densities array of densities at t +dt
	 * @param oldDensities array of densities at t
	 * @param uAra velocity array being applied in x dir
	 * @param vAra velocity array being applied in y dir
	 * 
	 * see also back and forth : bfecc error correction byungmoon kim
	 * see also vorticity confinement : ron fedkin
	 */	
	public void calcAdvectionJS(int bound, float[] densities, float[] oldDensities, float[] uAra, float[] vAra){
		int i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, deltaT0;
		deltaT0 = this.deltaT *  this.numCellsX;
		for (int i = 1 ; i <= this.numCellsX ; ++i) {
			for (int j = 1 ; j <= this.numCellsY ; ++j ) {
				x = i - deltaT0 * uAra[IX(i,j)]; 
				y = j - deltaT0 * vAra[IX(i,j)];
				//center of cells
				if (x < 0.5f) { x = 0.5f;} 
				if (x > this.numCellsX + 0.5) {	x = this.numCellsX + 0.5f;} 
				i0 = (int)x; 
				i1 = i0 + 1;
				
				if (y<0.5) {y=0.5f;} 
				if (y > this.numCellsY + 0.5){y = this.numCellsY + 0.5f;}
				j0=(int)y; j1=j0+1;
				
				s1 = x-i0; 
				s0 = 1-s1; 
				t1 = y-j0; 
				t0 = 1-t1;
				
				densities[IX(i,j)] = s0 * (t0 * oldDensities[IX(i0,j0)] + t1 * oldDensities[IX(i0,j1)]) 
								   + s1 * (t0 * oldDensities[IX(i1,j0)] + t1 * oldDensities[IX(i1,j1)]);
			}//for j
		}//for i
		this.setBoundariesJS(bound, densities);
	}//calcAdvection
	
	/**
	 * adds sources of density (mouseclicks, for example) to the density array
	 * @param dAra density ara
	 * @param sourceAra locations of mouseclicks
	 */
	public void addSource(float[] dAra, float[] sourceAra){ 
		for (int i=0 ; i < this.araSize ; i++ ) dAra[i] += this.deltaT * sourceAra[i]; 
		}
	
	/**
	 * calculates the appropriate pathways for the density through the velocity field
	 * @param xAra densities at time t+dt
	 * @param xOldAra old densities
	 * @param uAra velocity ara in x
	 * @param vAra velocity ara in y
	 */	
	public void calcDensityStep (float[] xAra, float[] xOldAra, float[] uAra, float[] vAra){
		this.addSource (xAra, xOldAra);
		int bound = 0;
		//SWAP ( x0, x )
			float[] tmpAra = xOldAra;
			xOldAra = xAra;
			xAra = tmpAra;
		this.calcDiffusionJS(bound, xAra, xOldAra, this.diffRate);
		
		//SWAP ( x0, x );
			tmpAra = xOldAra;
			xOldAra = xAra;
			xAra = tmpAra;
		this.calcAdvectionJS(bound, xAra, xOldAra, uAra, vAra);
	}//densityCalculationStep
	
	/**
	 * velocity solver
	 * @param uAra x dir velocity
	 * @param vAra y dir velocity
	 * @param uOldAra old velocity x
	 * @param vOldAra old velocity y
	 * @param visc viscosity
	 */
	public void calcVelocityStep(float[] uAra, float[] vAra, float[] uOldAra, float[] vOldAra){
		this.addSource( uAra, uOldAra); 
		this.addSource ( vAra, vOldAra);
//		float[] tmp = uOldAra;
//		uOldAra = uAra;
//		uAra = tmp;			
//	calcDiffusionJS(1, uAra, uOldAra, this.visc);
//		//SWAP ( vOldAra, vAra ); 
//		tmp = vOldAra;
//		vOldAra = vAra;
//		vAra = tmp;			
//	calcDiffusionJS (2, vAra, vOldAra, this.visc);
//	
//	this.project(uAra, vAra, uOldAra, vOldAra);
		
		//SWAP ( uOldAra, uAra ); 
			float[] tmp = uOldAra;
			uOldAra = uAra;
			uAra = tmp;			
		calcDiffusionJS(1, uAra, uOldAra, this.visc);
		
		//SWAP ( vOldAra, vAra ); 
			tmp = vOldAra;
			vOldAra = vAra;
			vAra = tmp;			
		calcDiffusionJS (2, vAra, vOldAra, this.visc);
		
		this.project(uAra, vAra, uOldAra, vOldAra);
		//SWAP ( uOldAra, uAra );
			tmp = uOldAra;
			uOldAra = uAra;
			uAra = tmp;			
		//SWAP ( vOldAra, vAra );
			tmp = vOldAra;
			vOldAra = vAra;
			vAra = tmp;			
	
		this.calcAdvectionJS(1, uAra, uOldAra, uOldAra, vOldAra); 
		this.calcAdvectionJS(2, vAra, vOldAra, uOldAra, vOldAra);
		//call project again to calculate pressures properly		
		this.project (uAra, vAra, uOldAra, vOldAra );
	}//velocityCalcStep
	
	/**
	 * cause a modification to the density array where clicked 
	 */
	 public void calcDensityEffect(int cellX, int cellY, float[] densityAra){
		 int xIDX, yIDX;
		 for(int i = -this.densityRadius; i <= this.densityRadius; ++i){
			 for(int j = -this.densityRadius; j <= this.densityRadius; ++j){
				 xIDX = i + cellX; 
				 yIDX = j + cellY;
				 xIDX += (xIDX < 0) ? this.numCellsX : 0;
				 xIDX -= (xIDX > this.numCellsX) ? this.numCellsX : 0;
				 yIDX += (yIDX < 0) ? this.numCellsY : 0;
				 yIDX -= (yIDX > this.numCellsY) ? this.numCellsY : 0;
				 densityAra[IX(xIDX, yIDX)] = this.densityStr * 10;
			 }//for j
		 }//for i
	 }//calcDensityEffect
	 
	/**
	 * process poisson calculation for preserving mass pressure projection
	 * @param uAra velocity field in x
	 * @param vAra velocity field in y
	 * @param pressure mass conserving field
	 * @param div gradient
	 * 
	 * better to use conjugate gradient
	 */	
	public void project(float[] uAra, float[] vAra, float[] pressure, float[] div){
		float height = 1.0f / this.numCellsX;
		for (int i = 1; i <= this.numCellsX; ++i ) {
			for (int j = 1; j <= this.numCellsY; ++j) {
				div[IX(i,j)] = -0.5f * height * (uAra[IX(i+1,j)]-uAra[IX(i-1,j)] + vAra[IX(i,j+1)] - vAra[IX(i,j-1)]);
				pressure[IX(i,j)] = 0;
			}//for j
		}//for i
		this.setBoundariesJS( 0, div ); 
		this.setBoundariesJS( 0, pressure );
		
		for (int k=0 ; k<20 ; k++ ) {
			for (int i = 1 ; i <= this.numCellsX ; ++i ) {
				for (int j = 1 ; j <= this.numCellsY ; ++j ) {
					pressure[IX(i,j)] = (div[IX(i,j)] + pressure[IX(i-1,j)] + pressure[IX(i+1,j)] + pressure[IX(i,j-1)] + pressure[IX(i,j+1)]) / 4.0f;
				}//for j
			}//for i
			this.setBoundariesJS( 0, pressure );
		}//for k
		
		for (int i = 1 ; i <= this.numCellsX ; ++i ) {
			for (int j = 1 ; j <= this.numCellsY ; ++j ) {
				uAra[IX(i,j)] -= 0.5f*(pressure[IX(i+1,j)] - pressure[IX(i-1,j)])/height;
				vAra[IX(i,j)] -= 0.5f*(pressure[IX(i,j+1)] - pressure[IX(i,j-1)])/height;
			}//for i
		}//for j
		this.setBoundariesJS(1, uAra); 
		this.setBoundariesJS(2, vAra);
	}//project

	/**
	 * stam solver index interpreter - takes x, y indices for 2 d array, returns index into 1 d array
	 * @param i x coord
	 * @param j y coord
	 * @return index into 1 d array
	 */
	public int IX(int i, int j){ return  i + ( this.numCellsX + 2 ) * j;}			

	//finds interpolation t of val between valMin and valMax and applies to mapMin and mapMax
	public float mapF(float val, float valMin, float valMax, float mapMin, float mapMax){
		if (Math.abs((valMax - valMin)) < this.epsVal) return 0;
		float minMaxDiff = valMax-valMin;
		float t = (val-valMin);		
		float calc = ((minMaxDiff-t) * mapMin + (t * mapMax)) / (1.0f * minMaxDiff);
		return calc;
	}//mapF

	/**
	 * print to console the current status of the program flags
	 */
	public void showFlags() {
		System.out.println();
		for (int i = 0; i < numFlags; ++i) {
			System.out.println(flagNames[i] + " = " + flags[i]);
		}
	}// showFlags function

	/**
	 * process keyboard input
	 */
	public void keyPressed(){
		switch (key){
			case 'i' :
			case 'I' : {//reinitialize sim
				this.initProgram();	
				showFlags();	
				break;
			}
			case ' ': {//start or stop simulation
				flags[runSim] = !flags[runSim];	
				showFlags();	
				break;
			}
			case 's' :
			case 'S' : {//save picture of current image
				save(sketchPath() + "\\" + ((flags[dispParticles]) ? ("partFluidImgs") : ("flowFluidImgs") )  + "\\" + "img" + (int)random(10) + "" + (int)random(10) + "" + (int) random(10) + "" + (int)random(10) + ((flags[dispParticles]) ? ("partFluid") : ("flowFluid") ) + ".jpg");
				break;
			}
			case 'a' :
			case 'A' : {//start/stop saving every frame for making into animation
				flags[saveAnim] = !flags[saveAnim];			
				break;
			}

			
			case 'c' : 
			case 'C' : {//toggle using custom or laplacian stencil for diffusion
				flags[dispParticles] = !flags[dispParticles];
				
				showFlags();	
				break;
			}

			case ',' :
			case '<' : {//decrease deltaT value to some lower bound
				deltaT -=.1f;
				deltaT = ((deltaT < 0) ? 0 : deltaT);
				System.out.println("Delta T = " +deltaT);
				break;
			}
			
			case '.' :
			case '>' : {//increase deltaT value to some upper bound
				deltaT +=.1f;
				deltaT = ((deltaT > 7) ? 7 : deltaT);
				System.out.println("Delta T = " +deltaT);
				break;
			}
	
			case ';' :
			case ':' : {//decrease the number of cycles between each draw, to some lower bound
				cycleModDraw -= 1;
				cycleModDraw = ((cycleModDraw < 1) ? 1 : cycleModDraw);
				System.out.println("Cycles to draw = " + cycleModDraw);
				break;
			}
	
			case '\'' :
			case '"' : {//increase the number of cycles between each draw to some upper bound
				cycleModDraw += 1;
				cycleModDraw = ((cycleModDraw > 20) ? 20 : cycleModDraw);
				System.out.println("Cycles to draw = " + cycleModDraw);
				break;
			}
			
			case '1' : {//set viscosity
				this.visc = 0;
				break;
			}
			
			case '2' :{
				this.visc = .0001f;
				break;
			}
	
			case '3' : {			
				this.visc = .001f;
				break;
			}
			
			case '4' : {
				this.visc = .01f;
				break;
			}		
			case '5' : {
				this.flags[concForce] = !this.flags[concForce]; //toggle concentric opposite forces
				break;
			}
							
			default : {				
			}
		}//switch	
	}//keypressed method
	public static void main(String[] passedArgs) {		
		String[] appletArgs = new String[] { "StamFluidSolverPkg.StamFluidSolverGlob" };
		if (passedArgs != null) {   	PApplet.main(PApplet.concat(appletArgs, passedArgs));    } else {   	PApplet.main(appletArgs);    }
	}//main


}// StamFluidSolverGlob class

