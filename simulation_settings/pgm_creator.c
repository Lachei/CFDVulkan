#include <stdio.h>
#include <stdlib.h>
  
    #define MAGIC_NUMBER "P2"
    #define MAX_VALUE 10
    #define TEST_CASE 0
    #define LEFT_BORDER 3
    #define RIGHT_BORDER 2
    #define TOP_BORDER 2
    #define BOTTOM_BORDER 0
    #define IMAGE_SIZE_X 2400
    #define IMAGE_SIZE_Y 800
    #define TREE_HEIGHT 0
    #define TREE_WIDTH 0
    #define TREE_POS 50
    #define TREES_DISTANCE 3
    #define TREE_TRUMK_HEIGHT 0 
    #define MOUNTAIN_HEIGHT 300
    #define MOUNTAIN_WIDTH_TOP  500
    #define MOUNTAIN_START_POS 100
    #define MOUNTAIN_SLOPE_STEPS_WIDTH_R 25
    #define MOUNTAIN_SLOPE_STEPS_HEIGHT_R 10
    #define MOUNTAIN_SLOPE_STEPS_WIDTH_L 10
    #define MOUNTAIN_SLOPE_STEPS_HEIGHT_L 20
    #define CHIMNEY_HEIGHT 10
    #define CHIMNEY_WIDTH 8
    #define HOUSES_START_POS 2000
    #define NUM_HOUSES 10
    #define HOUSE_SIZE 10
    #define HOUSE_DISTANCE 10

    void make_tree(int x_pos, int y_pos,int **image){
        for(int i=x_pos;i<x_pos+TREE_WIDTH;i++){
            for(int j=y_pos;j>y_pos-TREE_HEIGHT; j--){
                if(j>y_pos-TREE_TRUMK_HEIGHT+1)image[j][i]=0;
                else image[j][i]=10;
            }
        }
    }

   int main(){

    size_t i, j=0;
    int temp = 0;
  
    // Suppose the 2D Array to be converted to Image is as given below 
    int **image=malloc(IMAGE_SIZE_Y*sizeof(int*));
    for(int p=0;p<IMAGE_SIZE_Y;p++)image[p]=malloc(IMAGE_SIZE_X*sizeof(int));
    printf("Setting up grid with borders...\n");
    // Forming a white square inside a black square
    // Play here to form your own images!
    for (i = 0; i < IMAGE_SIZE_Y; ++i){
        for (j = 0; j < IMAGE_SIZE_X; ++j){
            if(j==0)image[i][j]=LEFT_BORDER;
            else if(j==IMAGE_SIZE_X-1)image[i][j]=RIGHT_BORDER;
            else if(i==0){
                if(j<(int)((1.0/5.0)*IMAGE_SIZE_X))image[i][j]=3;
                else image[i][j]=TOP_BORDER;
            }
            else if(i==IMAGE_SIZE_Y-1)image[i][j]=BOTTOM_BORDER;
            else image[i][j]=4;
        }
    }

    if(TEST_CASE==0){
        printf("Creating Mountain...\n");
        //MOUNTAIN
        size_t x_base = MOUNTAIN_START_POS;
        size_t y_base = IMAGE_SIZE_Y-1;
        j=x_base;
        size_t num_steps_l=0;
        size_t num_steps_r=0;
        printf("left slope...\n");
        while(num_steps_l<MOUNTAIN_HEIGHT/MOUNTAIN_SLOPE_STEPS_HEIGHT_L){
            for(i=y_base;i>y_base-num_steps_l*MOUNTAIN_SLOPE_STEPS_HEIGHT_L;i--){
                for(j=x_base+num_steps_l*MOUNTAIN_SLOPE_STEPS_WIDTH_L;j<x_base+num_steps_l*MOUNTAIN_SLOPE_STEPS_WIDTH_L+MOUNTAIN_SLOPE_STEPS_WIDTH_L;j++){
                    image[i][j]=0;
                }
            }
            num_steps_l++;
        }
        printf("mountain top...\n");
        size_t start_top=j;
        for(i=y_base;i>y_base-MOUNTAIN_HEIGHT;i--){
            for(j=start_top;j<start_top+MOUNTAIN_WIDTH_TOP;j++){
                image[i][j]=0;
                if(i==y_base-MOUNTAIN_HEIGHT+1 && j==start_top+MOUNTAIN_WIDTH_TOP/2){
                    printf("creating chimney...\n");
                    for(int y=y_base-MOUNTAIN_HEIGHT;y>y_base-MOUNTAIN_HEIGHT-CHIMNEY_HEIGHT;y--){
                        for(int x=j-CHIMNEY_WIDTH/2; x<j+CHIMNEY_WIDTH/2;x++){

                            if(y==y_base-MOUNTAIN_HEIGHT-CHIMNEY_HEIGHT+1 && x>j-CHIMNEY_WIDTH/2+1&&x<j+CHIMNEY_WIDTH/2-2){
                                image[y][x]=5;
                            }else image[y][x]=0;
                        }
                    }
                }
            }
        }
        printf("right slope with trees...\n");
        size_t start_slope_l=j;  
        while(num_steps_r<MOUNTAIN_HEIGHT/MOUNTAIN_SLOPE_STEPS_HEIGHT_R){
            for(i=y_base;i>y_base-MOUNTAIN_HEIGHT+num_steps_r*MOUNTAIN_SLOPE_STEPS_HEIGHT_R;i--){
                for(j=start_slope_l+num_steps_r*MOUNTAIN_SLOPE_STEPS_WIDTH_R;j<start_slope_l+num_steps_r*MOUNTAIN_SLOPE_STEPS_WIDTH_R+MOUNTAIN_SLOPE_STEPS_WIDTH_R;j++){
                    image[i][j]=0;
                    if(i==y_base-MOUNTAIN_HEIGHT+num_steps_r*MOUNTAIN_SLOPE_STEPS_HEIGHT_R+1 && j>start_slope_l+num_steps_r*MOUNTAIN_SLOPE_STEPS_WIDTH_R && (j-start_slope_l+num_steps_r*MOUNTAIN_SLOPE_STEPS_WIDTH_R)%(TREES_DISTANCE+TREE_WIDTH)==0){
                        make_tree(j,i,image);
                    }
                }
            }
            num_steps_r++;      
        }

        printf("creating houses...\n");
        for(int i=HOUSES_START_POS; i<HOUSES_START_POS+NUM_HOUSES*(HOUSE_DISTANCE+HOUSE_SIZE);i+=(HOUSE_DISTANCE+HOUSE_SIZE)){
            for(int x=i; x<i+HOUSE_SIZE; x++){
                for(int y=IMAGE_SIZE_Y-2; y>IMAGE_SIZE_Y-2-HOUSE_SIZE; y--){
                    image[y][x]=0;
                }
            }
        }
    }else if(TEST_CASE ==1){
        make_tree(TREE_POS, IMAGE_SIZE_Y-2, image);
    }

    image[0][0]=image[IMAGE_SIZE_Y-1][0]=image[IMAGE_SIZE_Y-1][IMAGE_SIZE_X-1]=image[0][IMAGE_SIZE_X-1]=0;   

    printf("writing to file\n");
    FILE* pgmimg; 

    //Creating a file ready to be written with a name of "myimg.pgm"
    pgmimg = fopen("death_valley_no_trees.pgm", "wb"); 
  
    // Writing Magic Number to the File 
    fprintf(pgmimg, "%s\n", MAGIC_NUMBER);  
  
    // Writing the size of the image 
    fprintf(pgmimg, "%d %d\n", IMAGE_SIZE_X, IMAGE_SIZE_Y);  
  
    // Writing the maximum gray value 
    fprintf(pgmimg, "%d\n", MAX_VALUE);

    int count = 0; 
    for (i = 0; i < IMAGE_SIZE_Y; i++) { 
        for (j = 0; j < IMAGE_SIZE_X; j++) { 
            temp = image[i][j]; 
  
            // Writing the gray values in the 2D array to the file 
            fprintf(pgmimg, "%d ", temp); 
        } 
        fprintf(pgmimg, "\n"); 
    } 
    fclose(pgmimg); 

    for(int p=0;p<IMAGE_SIZE_Y;p++)free(image[p]);
    free(image);
    printf("pgm creation finished. Success!\n");
} 