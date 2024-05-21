#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define SIZE 10
#define GRIDSIZE 60
#define BOUNDARY 100
#define pi 3.14159265358979

struct square{
    float angle;
    double xCenter;
    double yCenter;
} square1, square2;

typedef struct square Square;

void emptyGrid(char[GRIDSIZE][GRIDSIZE]);
void printGrid(char[GRIDSIZE][GRIDSIZE]);
void setCoordinate(char, char[GRIDSIZE][GRIDSIZE], double, double);
void setSquare(Square*, char[GRIDSIZE][GRIDSIZE]);
void emptySquare(Square*, char[GRIDSIZE][GRIDSIZE]);
float cot(float);
float rTan(float);
int orientation(float, float, float, float, float, float);
bool intersecting(float, float, float, float, float, float, float, float);
int separation(Square*, Square*, int*);
void cornerCoord(Square*, float[2], int);
void closestPoint(Square*, int*, float[2], float[2]);
float dot(float[2], float[2]);

int main(){
    char grid[GRIDSIZE][GRIDSIZE];
//fill grid with spaces
    emptyGrid(grid);
//initial values for s1 and s2

// sample coords 1
    // double s1positionX = 10;
    // double s1positionY = 10;
    // double s1velocityX = 0.075;
    // double s1velocityY = 0.075;
    // double s1anglePos = 0;
    // double s1angleVelocity = 0.23;

    // double s2positionX = 40;
    // double s2positionY = 35;
    // double s2velocityX = -0.075;
    // double s2velocityY = -0.075;
    // double s2anglePos = 67.5;
    // double s2angleVelocity = 0;

    // float e = -0.1; //elasticity
    // float m1 = 1; //mass of object 1
    // float m2 = 1; //mass of object 2
    // float moment1 = (m1*SIZE*SIZE)/3; //moment of inertia of object 1
    // float moment2 = (m2*SIZE*SIZE)/3; //moment of inertia of object 2

//sample coords 2

    double s1positionX = 8;
    double s1positionY = 35;
    double s1velocityX = 0.06;
    double s1velocityY = -0.02;
    double s1anglePos = 15;
    double s1angleVelocity = 0.05;

    double s2positionX = 58;
    double s2positionY = 25;
    double s2velocityX = -0.07;
    double s2velocityY = 0.03;
    double s2anglePos = 67.5;
    double s2angleVelocity = 0.022;

    float e = -0.8; //elasticity
    float m1 = 1; //mass of object 1
    float m2 = 2; //mass of object 2
    float moment1 = (m1*SIZE*SIZE)/3; //moment of inertia of object 1
    float moment2 = (m2*SIZE*SIZE)/3; //moment of inertia of object 2

//sample coords 3

    // double s1positionX = 3;
    // double s1positionY = 5;
    // double s1velocityX = 0.1;
    // double s1velocityY = 0.075;
    // double s1anglePos = 24;
    // double s1angleVelocity = 0.01;

    // double s2positionX = 58;
    // double s2positionY = 25;
    // double s2velocityX = -0.12;
    // double s2velocityY = 0;
    // double s2anglePos = 36;
    // double s2angleVelocity = 0.02;

    // float e = -0.5; //elasticity
    // float m1 = 1; //mass of object 1
    // float m2 = 2; //mass of object 2
    // float moment1 = (m1*SIZE*SIZE)/3; //moment of inertia of object 1
    // float moment2 = (m2*SIZE*SIZE)/3; //moment of inertia of object 2


//declaration of function things
    Square* s1 = &square1;
    Square* s2 = &square2;

    square1.angle = (pi/180)*(s1anglePos);
    square1.xCenter = s1positionX;
    square1.yCenter = s1positionY;
    square2.angle = (pi/180)*(s2anglePos);
    square2.xCenter = s2positionX;
    square2.yCenter = s2positionY;
//Rotate a square
    int frames = 400;
    int iFrames = 0;
    for(int i = 0; i < frames; i++){
        if(iFrames > 0){
            iFrames--;
        }
        // if(iFrames == 2){
        //     // printf("\n%f\n", s2angleVelocity);
        // }
        //Rotation and movement of Square 1
        square1.angle+=(pi/180)*s1angleVelocity;
        square1.xCenter+=s1velocityX;
        square1.yCenter+=s1velocityY;
        square2.angle+=(pi/180)*s2angleVelocity;
        square2.xCenter+=s2velocityX;
        square2.yCenter+=s2velocityY;
        if(square1.angle == 0){
            square1.angle = 0.00005;
        }
        if(square2.angle == 0){
            square2.angle = 0.00005;
        }
        //placement of squares
        setSquare(s1, grid);
        setSquare(s2,grid);
        printGrid(grid);
        // printf("(%f,%f)", s1->xCenter, s1->yCenter);
        // printf("(%f,%f)", s1velocityX, s1velocityY);
        // printf("%f & %f", s1anglePos, s1angleVelocity);
        int edge = 0;
        int *eP = &edge;
        int sep = separation(s1,s2,eP);
        if(sep > 0 && iFrames == 0){
            // printf("collision detected, stopping\n");
            // printf("Corner %d\n", sep);
            // printf("Frame #%d\n", i);
            iFrames = 10;
            float point[2];
            float closestP[2];
            if(sep/4 == 0){
                cornerCoord(s1, point, (sep%4));
                closestPoint(s2, eP, point, closestP);
                
            }else if(sep/4 == 1){
                cornerCoord(s2, point, (sep%4));
                closestPoint(s1, eP, point, closestP);
            }
            //printf("(%f, %f)", point[0], closestP[0]);
            //printf("(%f, %f)", point[1], closestP[1]);
            


            //Collision Behemoth & Evil Vector Magic
            //radius vectors
            // printf("\n\n%f\n\n",s2velocityY);
            float raP[2] = {closestP[0]-(s1->xCenter), closestP[1]-(s1->yCenter)};
            float rbP[2] = {closestP[0]-(s2->xCenter), closestP[1]-(s2->yCenter)};
            //initial velocity vectors at centers
            float vai[2] = {s2velocityX, s2velocityY};
            float vbi[2] = {s1velocityX, s1velocityY};
            //orthogonals to radius vectors (necessary for rotation)
            float raPortho[2] = {(raP[0]*cos(s1->angle)-raP[1]*sin(s1->angle)), (raP[0]*sin(s1->angle)+raP[1]*cos(s2->angle))};
            float rbPortho[2] = {(rbP[0]*cos(s2->angle)-rbP[1]*sin(s2->angle)), (rbP[0]*sin(s2->angle)+rbP[1]*cos(s2->angle))};
            //initial velocity vectors at points of contact
            float vapi[2] = {vai[0]+s1angleVelocity*raPortho[0], vai[1]+s1angleVelocity*raPortho[1]};
            float vbpi[2] = {vbi[0]+s2angleVelocity*rbPortho[0], vbi[1]+s2angleVelocity*rbPortho[1]};
            //total initial velocity vectors at points p
            float vabi[2] = {vapi[0]-vbpi[0],vapi[1]-vbpi[1]};
            //normal vector & copy
            float n1[2] = {point[0]-closestP[0], point[1]-closestP[1]};
            float n2[2] = {point[0]-closestP[0], point[1]-closestP[1]};
            //magic constant
            float k = -((1+e)*(dot(vabi,n1)))/((dot(n1,n2)*((1/m1)+(1/m2)))+((dot(raPortho,n1))*(dot(raPortho,n1)))/moment1+((dot(rbPortho,n1))*dot(rbPortho,n1))/moment2); //what the frick?
            //please be right
            s1angleVelocity+=((pi/180)*(dot(raPortho,n1))*(k/moment1));
            s2angleVelocity+=((pi/180)*(dot(rbPortho,n1))*(k/moment2));
            s1velocityX+=(n1[0]*(k/m1));
            s1velocityY+=(n1[1]*(k/m1));
            s2velocityX-=(n1[0]*(k/m2));
            s2velocityY-=(n1[1]*(k/m2));
            printf("<%f,%f>\n", n1[0], n1[1]);
            // printf("\n\n%f\n\n",s2velocityY);
            // printf("\n\n%f\n\n", k);
            // printf("%f,%f",);
            // break;
        }
        //frame reset
        usleep(4000);
        system("clear");
        emptySquare(s1, grid);
        emptySquare(s2,grid);
    }
    return 0;
}

void emptyGrid(char g[GRIDSIZE][GRIDSIZE]){
   for(int i = 0; i < GRIDSIZE; i++){
        for(int j = 0; j < GRIDSIZE; j++){
            g[i][j]='.';
        }
    }
}

void printGrid(char g[GRIDSIZE][GRIDSIZE]){
    for(int i = 0; i < GRIDSIZE; i++){
        int a = GRIDSIZE-1-i;
        for(int j = 0; j < GRIDSIZE; j++){
            printf("%c", g[a][j]);
        }
        printf("\n");
    }
}

void setCoordinate(char c, char g[GRIDSIZE][GRIDSIZE], double x, double y){
    int b = x;
    int a = y;
    if(!(g[a][b]=='0')){
        g[a][b]=c;
    }
}

void setSquare(Square* s, char g[GRIDSIZE][GRIDSIZE]){
    double xCorner = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle));
    double yCorner = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle));
    for(float i = 0; i < SIZE; i+=0.25){
        for(float j = 0; j < SIZE; j+=0.25){
            double x = xCorner+i*(cos(s->angle))-j*sin(s->angle);
            double y = yCorner+i*sin(s->angle)+j*cos(s->angle);
            if(x < 0 || y < 0 || x > GRIDSIZE - 1|| y > GRIDSIZE - 1){
            }else{
                setCoordinate('#', g, x, y);  
            }
            // setCoordinate('0', g, xCorner, yCorner);
        }
    }
}

void emptySquare(Square* s, char g[GRIDSIZE][GRIDSIZE]){
    double xCorner = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle));
    double yCorner = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle));
    for(float i = 0; i < SIZE; i+=0.25){
        for(float j = 0; j < SIZE; j+=0.25){
            double x = xCorner+i*(cos(s->angle))-j*sin(s->angle);
            double y = yCorner+i*sin(s->angle)+j*cos(s->angle);
            if(x < 0 || y < 0 || x > GRIDSIZE - 1|| y > GRIDSIZE - 1){
            }else{
                setCoordinate('.', g, x, y);  
            }
        }
    }
}

float cot(float theta){
    if(sin(theta) == 0){
        theta = 0.00005;
    }
    return (cos(theta)/sin(theta));
}

float rTan(float theta){
    if(cos(theta) == 0){
        theta = (pi/180)*90.00005;
    }
    return (sin(theta)/cos(theta));
}

int orientation(float ax, float ay, float bx, float by, float cx, float cy){
    float a = ((by-ay)*(cx-bx)-(bx-ax)*(cy-by)); //Written on geeksforgeeks
    if(a > 0){
        return 1;
    }else if(a < 0){
        return 2;
    }
    return 0;
};

bool intersecting(float ax, float ay, float bx, float by, float cx, float cy, float dx, float dy){
    int a = orientation(ax,ay,bx,by,cx,cy); //Written on geeksforgeeks
    int b = orientation(ax,ay,bx,by,dx,dy);
    int c = orientation(cx,cy,dx,dy,ax,ay);
    int d = orientation(cx,cy,dx,dy,bx,by);

    if(a != b && c != d){
        return true;
    }
    return false;
}

int separation(Square* s1, Square* s2, int* e){
    //shape 1 is key
    //condition 1 for shape 1
    double s1xCorner = (s1->xCenter)-(SIZE/2)*sin((pi/4)-(s1->angle));
    double s1yCorner = (s1->yCenter)-(SIZE/2)*cos((pi/4)-(s1->angle));
    double s2xCorner = (s2->xCenter)-(SIZE/2)*sin((pi/4)-(s2->angle));
    double s2yCorner = (s2->yCenter)-(SIZE/2)*cos((pi/4)-(s2->angle));
    //will check which corner of which square has collided
    int s1c1Check = 0;
    int s1c2Check = 0;
    int s1c3Check = 0;
    int s1c4Check = 0;
    int s2c1Check = 0;
    int s2c2Check = 0;
    int s2c3Check = 0;
    int s2c4Check = 0;


    float s1p1c1 = (s1xCorner)+cot(s1->angle)*(BOUNDARY-(s1yCorner));
    float s1p2c1 = (s1xCorner)-(SIZE*sin(s1->angle))+cot(s1->angle)*(BOUNDARY-(s1yCorner)-(SIZE*cos(s1->angle)));
    float s2p1c1 = (s2xCorner)+cot(s1->angle)*(BOUNDARY-(s2yCorner));
    float s2p2c1 = (s2xCorner)-(SIZE*sin(s2->angle))+cot(s1->angle)*(BOUNDARY-(s2yCorner)-(SIZE*cos(s2->angle)));
    float s2p3c1 = (s2xCorner)+(SIZE*cos(s2->angle))+cot(s1->angle)*(BOUNDARY-(s2yCorner)-(SIZE*sin(s2->angle)));
    float s2p4c1 = (s2xCorner)+(SIZE*(cos(s2->angle)-sin(s2->angle)))+cot(s1->angle)*(BOUNDARY-(s2yCorner)-(SIZE*(cos(s2->angle)+sin(s2->angle))));
    // printf("c1: %f\n%f\n%f\n%f\n%f\n%f\n\n\n", s1p1c1, s1p2c1, s2p1c1, s2p2c1, s2p3c1, s2p4c1);
    if(!(s1p1c1 < s2p1c1 && s2p1c1 < s1p2c1 ||
    s1p1c1 > s2p1c1 && s2p1c1 > s1p2c1 ||
    s1p1c1 < s2p2c1 && s2p2c1 < s1p2c1 ||
    s1p1c1 > s2p2c1 && s2p2c1 > s1p2c1 ||
    s1p1c1 < s2p3c1 && s2p3c1 < s1p2c1 ||
    s1p1c1 > s2p3c1 && s2p3c1 > s1p2c1 ||
    s1p1c1 < s2p4c1 && s2p4c1 < s1p2c1 ||
    s1p1c1 > s2p4c1 && s2p4c1 > s1p2c1)){
        return 0;
    }
    //Corner Check
    if(s1p1c1 < s2p1c1 && s2p1c1 < s1p2c1 ||
    s1p1c1 > s2p1c1 && s2p1c1 > s1p2c1){
        s2c1Check++;
    }
    if(s1p1c1 < s2p2c1 && s2p2c1 < s1p2c1 ||
    s1p1c1 > s2p2c1 && s2p2c1 > s1p2c1){
        s2c2Check++;
    }
    if(s1p1c1 < s2p3c1 && s2p3c1 < s1p2c1 ||
    s1p1c1 > s2p3c1 && s2p3c1 > s1p2c1){
        s2c3Check++;
    }
    if(s1p1c1 < s2p4c1 && s2p4c1 < s1p2c1 ||
    s1p1c1 > s2p4c1 && s2p4c1 > s1p2c1){
        s2c4Check++;
    }
    // printf("%d, %d, %d, %d, %d, %d, %d, %d\n\n\n", s1c1Check, s1c2Check, s1c3Check, s1c4Check, s2c1Check, s2c2Check, s2c3Check, s2c4Check);
    //condition 2 for shape 1
    float s1p1c2 = (s1xCorner)-rTan(s1->angle)*(BOUNDARY-(s1yCorner));
    float s1p2c2 = (s1xCorner)+(SIZE*cos(s1->angle))-rTan(s1->angle)*(BOUNDARY-(s1yCorner)-SIZE*sin(s1->angle));
    float s2p1c2 = (s2xCorner)-rTan(s1->angle)*(BOUNDARY-(s2yCorner));
    float s2p2c2 = (s2xCorner)-(SIZE*sin(s2->angle))-rTan(s1->angle)*(BOUNDARY-(s2yCorner)-SIZE*cos(s2->angle));
    float s2p3c2 = (s2xCorner)+(SIZE*cos(s2->angle))-rTan(s1->angle)*(BOUNDARY-(s2yCorner)-SIZE*sin(s2->angle));
    float s2p4c2 = (s2xCorner)-(SIZE*(sin(s2->angle)-cos(s2->angle)))-rTan(s1->angle)*(BOUNDARY-(s2yCorner)-SIZE*(sin(s2->angle)+cos(s2->angle)));
    // printf("c2:%f\n%f\n%f\n%f\n%f\n%f\n\n\n", s1p1c2, s1p2c2, s2p1c2, s2p2c2, s2p3c2, s2p4c2);
    if(!(s1p1c2 < s2p1c2 && s2p1c2 < s1p2c2 ||
    s1p1c2 > s2p1c2 && s2p1c2 > s1p2c2 ||
    s1p1c2 < s2p2c2 && s2p2c2 < s1p2c2 ||
    s1p1c2 > s2p2c2 && s2p2c2 > s1p2c2 ||
    s1p1c2 < s2p3c2 && s2p3c2 < s1p2c2 ||
    s1p1c2 > s2p3c2 && s2p3c2 > s1p2c2 ||
    s1p1c2 < s2p4c2 && s2p4c2 < s1p2c2 ||
    s1p1c2 > s2p4c2 && s2p4c2 > s1p2c2)){
        return 0;
    }
    //Corner Check
    if(s1p1c2 < s2p1c2 && s2p1c2 < s1p2c2 ||
    s1p1c2 > s2p1c2 && s2p1c2 > s1p2c2 ){
        s2c1Check++;
    }
    if(s1p1c2 < s2p2c2 && s2p2c2 < s1p2c2 ||
    s1p1c2 > s2p2c2 && s2p2c2 > s1p2c2){
        s2c2Check++;
    }
    if(s1p1c2 < s2p3c2 && s2p3c2 < s1p2c2 ||
    s1p1c2 > s2p3c2 && s2p3c2 > s1p2c2){
        s2c3Check++;
    }
    if(s1p1c2 < s2p4c2 && s2p4c2 < s1p2c2 ||
    s1p1c2 > s2p4c2 && s2p4c2 > s1p2c2){
        s2c4Check++;
    }
    // printf("%d, %d, %d, %d, %d, %d, %d, %d\n\n\n", s1c1Check, s1c2Check, s1c3Check, s1c4Check, s2c1Check, s2c2Check, s2c3Check, s2c4Check);
    //shape 2 is now key
    //condition 3 for shape 2
    float s2p1c3 = (s2xCorner)+cot(s2->angle)*(BOUNDARY-(s2yCorner));
    float s2p2c3 = (s2xCorner)-(SIZE*sin(s2->angle))+cot(s2->angle)*(BOUNDARY-(s2yCorner)-(SIZE*cos(s2->angle)));
    float s1p1c3 = (s1xCorner)+cot(s2->angle)*(BOUNDARY-(s1yCorner));
    float s1p2c3 = (s1xCorner)-(SIZE*sin(s1->angle))+cot(s2->angle)*(BOUNDARY-(s1yCorner)-(SIZE*cos(s1->angle)));
    float s1p3c3 = (s1xCorner)+(SIZE*cos(s1->angle))+cot(s2->angle)*(BOUNDARY-(s1yCorner)-(SIZE*sin(s1->angle)));
    float s1p4c3 = (s1xCorner)+(SIZE*(cos(s1->angle)-sin(s1->angle)))+cot(s2->angle)*(BOUNDARY-(s1yCorner)-(SIZE*(cos(s1->angle)+sin(s1->angle))));
    // printf("%f\n%f\n%f\n%f\n%f\n%f\n\n\n", s2p1c3, s2p2c3, s1p1c3, s1p2c3, s1p3c3, s1p4c3);
    if(!(s2p1c3 < s1p1c3 && s1p1c3 < s2p2c3 ||
    s2p1c3 > s1p1c3 && s1p1c3 > s2p2c3 ||
    s2p1c3 < s1p2c3 && s1p2c3 < s2p2c3 ||
    s2p1c3 > s1p2c3 && s1p2c3 > s2p2c3 ||
    s2p1c3 < s1p3c3 && s1p3c3 < s2p2c3 ||
    s2p1c3 > s1p3c3 && s1p3c3 > s2p2c3 ||
    s2p1c3 < s1p4c3 && s1p4c3 < s2p2c3 ||
    s2p1c3 > s1p4c3 && s1p4c3 > s2p2c3)){
        return 0;
    }
    //Corner Check
    if(s2p1c3 < s1p1c3 && s1p1c3 < s2p2c3 ||
    s2p1c3 > s1p1c3 && s1p1c3 > s2p2c3){
        s1c1Check++;
    }
    if(s2p1c3 < s1p2c3 && s1p2c3 < s2p2c3 ||
    s2p1c3 > s1p2c3 && s1p2c3 > s2p2c3){
        s1c2Check++;
    }
    if(s2p1c3 < s1p3c3 && s1p3c3 < s2p2c3 ||
    s2p1c3 > s1p3c3 && s1p3c3 > s2p2c3){
        s1c3Check++;
    }
    if(s2p1c3 < s1p4c3 && s1p4c3 < s2p2c3 ||
    s2p1c3 > s1p4c3 && s1p4c3 > s2p2c3){
        s1c4Check++;
    }
    // printf("%d, %d, %d, %d, %d, %d, %d, %d\n\n\n", s1c1Check, s1c2Check, s1c3Check, s1c4Check, s2c1Check, s2c2Check, s2c3Check, s2c4Check);
    //condition 4 for shape 2
    float s2p1c4 = (s2xCorner)-rTan(s2->angle)*(BOUNDARY-(s2yCorner));
    float s2p2c4 = (s2xCorner)+(SIZE*cos(s2->angle))-rTan(s2->angle)*(BOUNDARY-(s2yCorner)-SIZE*sin(s2->angle));
    float s1p1c4 = (s1xCorner)-rTan(s2->angle)*(BOUNDARY-(s1yCorner));
    float s1p2c4 = (s1xCorner)-(SIZE*sin(s1->angle))-rTan(s2->angle)*(BOUNDARY-(s1yCorner)-SIZE*cos(s1->angle));
    float s1p3c4 = (s1xCorner)+(SIZE*cos(s1->angle))-rTan(s2->angle)*(BOUNDARY-(s1yCorner)-SIZE*sin(s1->angle));
    float s1p4c4 = (s1xCorner)-(SIZE*(sin(s1->angle)-cos(s1->angle)))-rTan(s2->angle)*(BOUNDARY-(s1yCorner)-SIZE*(sin(s1->angle)+cos(s1->angle)));
    // printf("%f\n%f\n%f\n%f\n%f\n%f\n\n\n", s2p1c4, s2p2c4, s1p1c4, s1p2c4, s1p3c4, s1p4c4);
    if(!(s2p1c4 < s1p1c4 && s1p1c4 < s2p2c4 ||
    s2p1c4 > s1p1c4 && s1p1c4 > s2p2c4 ||
    s2p1c4 < s1p2c4 && s1p2c4 < s2p2c4 ||
    s2p1c4 > s1p2c4 && s1p2c4 > s2p2c4 ||
    s2p1c4 < s1p3c4 && s1p3c4 < s2p2c4 ||
    s2p1c4 > s1p3c4 && s1p3c4 > s2p2c4 ||
    s2p1c4 < s1p4c4 && s1p4c4 < s2p2c4 ||
    s2p1c4 > s1p4c4 && s1p4c4 > s2p2c4)){
        return 0;
    }
    //Corner Check
    if(s2p1c4 < s1p1c4 && s1p1c4 < s2p2c4 ||
    s2p1c4 > s1p1c4 && s1p1c4 > s2p2c4){
        s1c1Check++;
    }
    if(s2p1c4 < s1p2c4 && s1p2c4 < s2p2c4 ||
    s2p1c4 > s1p2c4 && s1p2c4 > s2p2c4){
        s1c2Check++;
    }
    if(s2p1c4 < s1p3c4 && s1p3c4 < s2p2c4 ||
    s2p1c4 > s1p3c4 && s1p3c4 > s2p2c4){
        s1c3Check++;
    }
    if(s2p1c4 < s1p4c4 && s1p4c4 < s2p2c4 ||
    s2p1c4 > s1p4c4 && s1p4c4 > s2p2c4){
        s1c4Check++;
    }
    // printf("%d, %d, %d, %d, %d, %d, %d, %d\n\n\n", s1c1Check, s1c2Check, s1c3Check, s1c4Check, s2c1Check, s2c2Check, s2c3Check, s2c4Check);
    //Edge finder
    float point1[2];
    float point2[2];
    cornerCoord(s1,point1,1);
    cornerCoord(s1,point2,2);
    if(s2c1Check==2 || s2c2Check==2||s2c3Check==2||s2c4Check==2){
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 1;
            //printf("Left edge, square 1\n");
        }
        cornerCoord(s1,point1,2);
        cornerCoord(s1,point2,3);
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 2;
           // printf("(%f, %f)\n(%f, %f)\n(%f, %f)\n(%f,%f)\n",s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1]);
           // printf("Top edge, square 1\n");       
        }
        cornerCoord(s1,point1,3);
        cornerCoord(s1,point2,0); 
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 3;
            //printf("Right edge, square 1\n");
        }
        cornerCoord(s1,point1,0);
        cornerCoord(s1,point2,1);
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 4;
            //printf("(%f, %f)\n(%f, %f)\n(%f, %f)\n(%f,%f)\n",s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1]);
            // printf("Bottom edge, square 1\n");
        }
    }else if(s1c1Check==2||s1c2Check==2||s1c3Check==2||s1c4Check==2){
        cornerCoord(s2,point1,1);
        cornerCoord(s2,point2,2);
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 5;
            // printf("Left edge, square 2\n");
        }
        cornerCoord(s2,point1,2);
        cornerCoord(s2,point2,3);
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 6;
            // printf("Top edge, square 2\n");       
        }
        cornerCoord(s2,point1,3);
        cornerCoord(s2,point2,0);         
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 7;
            // printf("Right edge, square 2\n");
        }
        cornerCoord(s2,point1,0);
        cornerCoord(s2,point2,1);
        if(intersecting(s1->xCenter,s1->yCenter,s2->xCenter,s2->yCenter,point1[0],point1[1],point2[0],point2[1])==true){
            *e = 8;
            // printf("Bottom edge, square 2\n");
        }
    }
    //Corner finder
    if(s1c1Check == 2){
        // printf("Bottom left, square 1\n");
        return 1;
    }else if(s1c2Check == 2){
        // printf("Top left, square 1\n");
        return 2;
    }else if(s1c3Check == 2){
        // printf("Bottom right, square 1\n");
        return 3;
    }else if(s1c4Check == 2){
        // printf("Top right, square 1\n");
        return 4;
    }else if(s2c1Check == 2){
        // printf("Bottom left, square 2\n");
        return 5;
    }else if(s2c2Check == 2){
        // printf("Top left, square 2\n");
        return 6;
    }else if(s2c2Check == 2){
        // printf("Bottom right, square 2\n");
        return 7;
    }else if(s2c2Check == 2){
        // printf("Top right, square 2\n");
        return 8;
    }
    return -1;
}

void cornerCoord(Square* s, float p[2], int c){
    if(c == 1){
        p[0] = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle));
        p[1] = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle));
    }else if(c == 2){
        p[0] = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle))-(SIZE*(sin(s->angle)));
        p[1] = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle))+(SIZE*(cos(s->angle)));
    }else if(c==3){
        p[0] = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle))+SIZE*(cos(s->angle));
        p[1] = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle))+SIZE*(sin(s->angle));
    }else if(c==0){
        p[0] = (s->xCenter)-(SIZE/2)*sin((pi/4)-(s->angle))-SIZE*(sin(s->angle)-cos(s->angle));
        p[1] = (s->yCenter)-(SIZE/2)*cos((pi/4)-(s->angle))+SIZE*(cos(s->angle)+sin(s->angle));
    }
}

void closestPoint(Square* s, int* e, float p[2], float p2[2]){
    float corner[2];
    cornerCoord(s, corner, ((*e)%4));
    // printf("Point and angle: (%f, %f), %f\n\n\n\n", corner[0], corner[1], (s->angle)*(180/pi));
    p2[0] = (((p[0])*(sin(s->angle))*(sin(s->angle)))+((corner[0])*cos(s->angle)*cos(s->angle))+(corner[1]*sin(s->angle)*cos(s->angle))-(p[1]*sin(s->angle)*cos(s->angle)));
    //((corner[0]-p[0]+(p[1]*rTan(s->angle))+(p[1]*cot(s->angle)))/(cot(s->angle)+rTan(s->angle)));
    p2[1] = ((p[1]*cos(s->angle)*cos(s->angle))-p[0]*cos(s->angle)*sin(s->angle)+corner[0]*sin(s->angle)*cos(s->angle)+corner[1]*sin(s->angle)*sin(s->angle));
    //((corner[1]-p[1]+(corner[0]*cot(s->angle))+(p[0]*rTan(s->angle)))/(cot(s->angle)+rTan(s->angle)));
}

float dot(float a[2], float b[2]){
    return (a[0]*b[0]+a[1]+b[1]);
}