#include<iostream>
#include<SFML/Graphics.hpp>
#include<vector>
#include<math.h>
#include<string.h> //use memset function

using namespace std;
using namespace sf;

//window size:
Vector2i wSize(1200, 1200);

//define mesh:
const int meshSize = 220, meshSizeP2 = 222; // square mesh, size excludes the boundary cells; meshSize is no. of cells in each direction
float cellWidth, cellHeight;

//Velocity scale and time step:
float Uscale = 200.f;
float dt = 0.05;

int cellIndxMap(int i, int j)   //i for vertical, j for horizontal
{
    // return the 1D equivalent index of cell, where each cell has 4 points
    return 4*(j + meshSize*i);
}

float lerp(float a, float b, float k)
{
    return a+k*(b-a);
}

void setBnd(int b, float field[meshSizeP2][meshSizeP2])
{
    for (int i = 1; i <= meshSize; i++)
    {
        field[0][i] = b==2 ? -field[1][i] : field[1][i];
        field[meshSize+1][i] = b==2 ? -field[meshSize][i] : field[meshSize][i];
        field[i][0] = b==1 ? -field[i][1] : field[i][1];
        field[i][meshSize+1] = b==1 ? -field[i][meshSize] : field[i][meshSize];
    
        // field[0][i] = b!=0 ? -field[1][i] : field[1][i];
        // field[meshSize+1][i] = b!=0 ? -field[meshSize][i] : field[meshSize][i];
        // field[i][0] = b!=0 ? -field[i][1] : field[i][1];
        // field[i][meshSize+1] = b!=0 ? -field[i][meshSize] : field[i][meshSize];
    
    }
    //set edges:
    field[0][0] = 0.5*(field[1][0] + field[0][1]);
    field[0][meshSize+1] = 0.5*(field[1][meshSize+1] + field[0][meshSize]);
    field[meshSize+1][0] = 0.5*(field[meshSize][0] + field[1][meshSize+1]);
    field[meshSize+1][meshSize+1] = 0.5*(field[meshSize][meshSize+1] + field[meshSize+1][meshSize]);
}

void diffuseField(int b, float field[meshSizeP2][meshSizeP2], float delta)   //pass 2d array by reference
{
    float fieldNext[meshSizeP2][meshSizeP2] = {};
    float del = delta;    //del = 4*kappa*dt/h^2, it controls the diffusion strength
    int iterMax = 2;   // should compare with prev iteration step to terminate iteration loop
    for (int iter = 0; iter <= iterMax; iter++)
    {
       for (int i = 0; i < meshSize; i++)
       {
           for (int j = 0; j < meshSize; j++)
           {
                if(iter==0) fieldNext[i+1][j+1] = field[i+1][j+1];  //start with prev field, boundary cells are already always zero
                // else if(iter<iterMax && iter>0)    fieldNext[i+1][j+1] = (fieldNext[i+1][j+1] + del*(1.0/4.0)*(field[i+1+1][j+1]+field[i+1-1][j+1]+field[i+1][j+1+1]+field[i+1][j+1-1] ))/(1+del);
                else if(iter<iterMax && iter>0)    fieldNext[i+1][j+1] = (field[i+1][j+1] + del*(1.0/4.0)*(fieldNext[i+1+1][j+1]+fieldNext[i+1-1][j+1]+fieldNext[i+1][j+1+1]+fieldNext[i+1][j+1-1] ))/(1+del);
                
              //  else if(iter<iterMax && iter>0)    fieldNext[i+1][j+1] = (field[i+1][j+1] + del*(1.0/4.0)*(field[i+1+1][j+1]+field[i+1-1][j+1]+field[i+1][j+1+1]+field[i+1][j+1-1] ))/(1+del);
                
                else field[i+1][j+1] = fieldNext[i+1][j+1];
           }
       }
      // if(b!=-1) setBnd(b, fieldNext);
    }
    if(b!=-1) setBnd(b, field);
    //memcpy(field, fieldNext, sizeof(field));
}

void advectField(int b, float field[][meshSizeP2], float field0[][meshSizeP2], float uxField[][meshSizeP2], float uyField[][meshSizeP2])
{
    float xPos, yPos, value1, value2;
    float tempField[meshSizeP2][meshSizeP2] = {};
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            xPos = (j)*cellWidth - uxField[i+1][j+1]*dt;
            yPos = (i)*cellHeight - uyField[i+1][j+1]*dt;
            int ix = floor(xPos/(float)cellWidth), iy = floor(yPos/(float)cellHeight);  // integer part
            float jx = (xPos/(float)cellWidth) - ix, jy = (yPos/(float)cellHeight) - iy;    //fractional part
           // ix = jx<=0.5?ix-1:ix;
           // iy = jy<=0.5?iy-1:iy;   //no need to change jx and jy
            if(ix+1<0 || ix+1>meshSize)  { cout<<"Problem in advection!"<<endl; ix=0; }   
            if(iy+1<0 || iy+1>meshSize)    iy=0;   //use zero values at boundaries 
            
            value1 = lerp(field0[iy+1][ix+1], field0[iy+1][ix+1+1], jx);  // second index for horizontal and first for vertical, be careful...
            value2 = lerp(field0[iy+1+1][ix+1], field0[iy+1+1][ix+1+1], jx);
            tempField[i+1][j+1] = lerp(value1, value2, jy);
        }
    }
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            field[i+1][j+1] = tempField[i+1][j+1];
        }
    }

    if(b!= -1)    setBnd(b, field);
    //memcpy(field, tempField, sizeof(field));
}

void project(float pField[meshSizeP2][meshSizeP2], float uxField[][meshSizeP2], float uyField[][meshSizeP2])
{
    int iterMax = 15;   // should compare with prev iteration step to terminate iteration loop
    float h = cellWidth; //same as cellHeight
    float divField[meshSizeP2][meshSizeP2] = {};
    for (int iter = 0; iter <= iterMax; iter++)
    {
       for (int i = 0; i < meshSize; i++)
       {
            for (int j = 0; j < meshSize; j++)
            {
                if(iter==0)    divField[i+1][j+1] = h*0.25*(uxField[i+1][j+1+1] - uxField[i+1][j+1-1] + uyField[i+1+1][j+1] - uyField[i+1-1][j+1] )/2.f;    // pField includes delT factor
                pField[i+1][j+1] = -1.0*divField[i+1][j+1] + 0.25*(pField[i+1+1][j+1]+pField[i+1-1][j+1]+pField[i+1][j+1+1]+pField[i+1][j+1-1] );
            }
       }
       //continuity for p at boundaries:
       setBnd(0, pField);
    }
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            uxField[i+1][j+1] -= (pField[i+1][j+1+1]-pField[i+1][j+1-1])/(2.f*h);
            uyField[i+1][j+1] -= (pField[i+1+1][j+1]-pField[i+1-1][j+1])/(2.f*h);
        }
    }
    setBnd(1, uxField);
    setBnd(2, uyField);
}

int main()
{
    //define fields:
    float density[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float ux[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uy[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uxPrev[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uyPrev[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    
    float pField[meshSizeP2][meshSizeP2] = {};

    cellWidth = (float)wSize.x/(float)meshSize;
    cellHeight = (float)wSize.y/(float)meshSize;

    dt = 2.f*cellWidth/Uscale;
    dt = dt>0.08?dt:0.08;
    cout<<"dt: "<<dt<<endl;
    

    RenderWindow w(VideoMode(wSize.x, wSize.y), "simplest2Dsolver");
    w.setFramerateLimit(60);
    //**********Create Mesh***********
    //define vertex array (quads) for each cell:
    vector<Vertex> cells; //4 vertices needed to define quad, need to be in a particular order (clockwise or anti-clockwise)  
    for (int i = 0; i < meshSize; i++)  //cellsPrev stores previous color info
    {
        for (int j = 0; j < meshSize; j++)
        {
            Vertex temp;
            temp.color = Color(255*density[i+1][j+1]);
            temp.position = Vector2f(j*cellWidth, i*cellHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f((j+1)*cellWidth, i*cellHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f((j+1)*cellWidth, (i+1)*cellHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f(j*cellWidth, (i+1)*cellHeight);
            cells.push_back(temp);
        }   
    }
    
    cout<<"Total no. of cells by meshSize squared: "<<cells.size()/(4*meshSize*meshSize)<<endl;   //should be 1
    cout<<"cellWidth: "<<cellWidth<<endl;
    //Use mouse click to add density
    Vector2i localMousePos(0, 0), localMousePrevPos(0, 0);
    Color currColor = Color::White;
    
    bool begin = false, mouseDown = false;
    float dissipation = 0.04;

    int ct = 0;
    int Niter=0;
    while (w.isOpen())
    {
        Event e;
        while (w.pollEvent(e))
        {
            if(e.type==Event::Closed)   w.close();
            if(e.type==Event::KeyPressed)   
            {
                if(e.key.code==Keyboard::C) 
                {
                    for (int i = 0; i < meshSize; i++)
                    {
                        for (int j = 0; j < meshSize; j++)
                        {
                            density[i+1][j+1] = 0.f;
                            ux[i+1][j+1] = 0.f;
                            uy[i+1][j+1] = 0.f;
                            uxPrev[i+1][j+1] = 0.f;
                            uyPrev[i+1][j+1] = 0.f;
                        }
                    }
                }
                if(e.key.code==Keyboard::B) begin=begin?false:true;
                if(e.key.code==Keyboard::E) Niter+=1;
                
                if(e.key.code==Keyboard::T) 
                {
                    int ix = (floor)((float)localMousePos.x/cellWidth);
                    int iy = (floor)((float)localMousePos.y/cellHeight);
                    float wz = (uyPrev[iy+1][ix+1+1] - uyPrev[iy+1][ix+1-1] - uxPrev[iy+1+1][ix+1] + uxPrev[iy+1-1][ix+1])/(2.f*cellWidth);
                    cout<<"vorticity:"<<wz<<endl;
                }

            }
            if(e.type==Event::MouseButtonReleased)
            {
                if(e.mouseButton.button==Mouse::Left)    mouseDown = false;
            }
            if(Mouse::isButtonPressed(Mouse::Left) && !mouseDown) 
            {   
                mouseDown = true;   
                localMousePrevPos = Mouse::getPosition(w); 
               // currColor = Color(255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX);
            }
        }

        //add source:
        localMousePos = Mouse::getPosition(w);
        int sourceSize = 5; //radius of source 
        if(mouseDown)
        {
            int localJ = floor(localMousePos.x/cellWidth), localI = floor(localMousePos.y/cellHeight);
            for (int iSource = localI-sourceSize; iSource <= localI+sourceSize; iSource++)
            {
                if(iSource<0 || iSource>=meshSize)   continue;
                for (int jSource = localJ-sourceSize; jSource <= localJ+sourceSize; jSource++)
                {
                    if(jSource<0 || jSource>=meshSize)   continue;
                    
                    uyPrev[iSource+1][jSource+1] = -Uscale;
                    // if(localMousePos.x != localMousePrevPos.x && localMousePos.y != localMousePrevPos.y)
                    // {
                    //     ux[iSource+1][jSource+1] = Uscale*tanh(localMousePos.x-localMousePrevPos.x);
                    //     uy[iSource+1][jSource+1] = Uscale*tanh(localMousePos.y-localMousePrevPos.y);
                    // }
                    // else    uy[iSource+1][jSource+1] = -Uscale;
                    float localD = abs(jSource - localJ)*cellWidth + abs(iSource - localI)*cellHeight, dMax = sourceSize*(cellWidth+cellHeight);
                    
                    density[iSource + 1][jSource + 1] = density[iSource + 1][jSource + 1]<1.f-localD/(2.f*dMax)?1.f-localD/(2.f*dMax):density[iSource + 1][jSource + 1]; // added ones for the boundary offset
                    int cellCoord = cellIndxMap(iSource, jSource);
                    for(int cIndx = 0; cIndx < 4; cIndx++)  cells[cellCoord + cIndx].color = currColor;
                }
                
            }
        }   

        if(begin)
        {
            
            //*******Calculate velocity for next instant:**********
           // diffuseField(ux, 0.001); // second argument specifies viscosity
           // diffuseField(uy, 0.001);

            // diffuseField(1, uxPrev, 1);
            // diffuseField(2, uyPrev, 1);

            // project(pField, uxPrev, uyPrev);
            
           // if(ct<=Niter)
            {
                advectField(1, ux, uxPrev, uxPrev, uyPrev);
                advectField(2, uy, uyPrev, uxPrev, uyPrev);
                
                project(pField, ux, uy);

                //diffuse:
                diffuseField(-1, density, 0.1); //second argument characterizes diffusion constant

                //advect:
                advectField(-1, density, density, ux, uy);
                ct++;
            }  
            for (int i = 0; i < meshSizeP2; i++)
            {
                for (int j = 0; j < meshSizeP2; j++)
                {
                    uxPrev[i][j] = ux[i][j];
                    uyPrev[i][j] = uy[i][j];
                    density[i][j] -= dissipation*dt*density[i][j]; //dissipate density
                }
            }

           // diffuseField(ux, 0.01); // second argument specifies viscosity
           // diffuseField(uy, 0.01);

//            project(pField, ux, uy);
            // //*****************
            //update density for cells:
            
            for (int i = 0; i < meshSize; i++)
            {
                for (int j = 0; j < meshSize; j++)
                {
                    int cellCoord = cellIndxMap(i, j);
                    
                    float tempVal = 0.25*(density[i][j]+density[i+1][j]+density[i][j+1]+density[i+1][j+1]);
                    //float tempVal = density[IX(i, j)]; 
                    Color tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord].color = tmpColor;
                    
                    tempVal = 0.25*(density[i][j+1]+density[i+1][j+1]+density[i][j+2]+density[i+1][j+2]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+1].color = tmpColor;

                    tempVal = 0.25*(density[i+1][j+1]+density[i+2][j+1]+density[i+1][j+2]+density[i+2][j+2]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+2].color = tmpColor;

                    tempVal = 0.25*(density[i+1][j]+density[i+2][j]+density[i+1][j+1]+density[i+2][j+1]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+3].color = tmpColor;
                }
            }
        
        }

                
        w.clear();
        w.draw(cells.data(), cells.size(), Quads);
        w.display();    
    
    
    }
    
    return 0;
}