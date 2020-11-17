# include <stdio.h>
# include <conio.h>
# include <math.h>
# include <stdlib.h>
#define Max 100
#define eps 0.0001
//Ham hoan doi
swap(float &a,float &b){
	float temp;
	temp=a;
	a=b;
	b=temp;
}
//xuat phuong trinh
void xuathe(float A[][Max],int n){
	printf("He phuong trinh can giai:\n");
	for(int i=0;i<n;i++){
		printf("%.3fX1 ",A[i][0]);
	    for(int j=1;j<n;j++){
	    	if(A[i][j]>=0) printf(" +");
	    	printf(" %.3fX%d",A[i][j],j+1);
		}
		printf("=");
		printf(" %.3f\n",A[i][n]);     
    }
}
// kiem tra he so aii co bang 0 khong neu co thi doi dong
void aiikhackhong(float A[][Max],int n){
	int i,j,k;
	for(i=0;i<n;i++){
		if(A[i][i]==0)
		   for(j=0;j<n;j++)
			   if(A[j][i]!=0&&A[i][j]!=0)
			       for(k=0;k<=n;k++)
			           swap(A[i][k],A[j][k]);
	    if(A[i][i]==0){
	    	printf("\n khong the tinh vi khong dua ve dang 2 !");
	        exit(1);
		}			   				    
	}
}
//Giai phuong trinh bang phuong phap giam du
int giamdu(float B[][Max],float R[],float X[],int n){
	printf("PHUONG PHAP GIAM DU:\n");
	int i,j,k,t;
	float temp,max,A[Max][Max],A2[Max][Max],R2[Max],X2[Max];
	for(i=0;i<n;i++){
	    for(j=0;j<=n;j++) A[i][j]=B[i][j];
	}
	//ham kiem tra aii co bang khong hay khong neu co thi doi dong
	aiikhackhong(A,n);
    printf("\nNhap vao vecto nghiem ban dau:\n  ");
    for(i=0;i<n;i++){
    	printf("x%d= ",i+1);
    	scanf("%f",&X[i]);
    	X2[i]=X[i];
	}
    // bien doi ve dang dang (2)
    for(i=0;i<n;i++){
	    temp=A[i][i];
        for(j=0;j<=n;j++){
        	A[i][j]/=temp;
        	A2[i][j]=A[i][j];
        }
    }     	
    
    // kiem tra co tinh duoc bang phuong phap giam du khong
    // tinh r[i] ban dau
    for(i=0;i<n;i++){
	    R[i]=A[i][n];
        for(j=0;j<n;j++){
             R[i]-=A[i][j]*X[j]	;
		}
		R2[i]=R[i];
    }
    int kt=0; // bien de kiem tra
	// Tim | rs | = max {|r[i]|} (i = 1->n) & tinh lai xs
	do{
		kt++;
		// tat vong lap khi lap qua nhieu
		if(kt>30){
			printf("\nkhong the tinh bang phuong phap nay!");
			return(0);
		} 
		t=0; // de thoat vong lap
		// tim max(|R[i]|)
		max=fabs(R[0]);
		k=0;
        for(i=0;i<n;i++){
        	if(max<fabs(R[i])){
			   max=fabs(R[i]);
        	   k=i;
            }
		}	
	    X[k]+=R[k];
		// tinh lai R[i] va kiem tra truong hop tiep theo
		temp=R[k];
		for(i=0;i<n;i++){
			R[i]-=A[i][k]*temp;
			if(fabs(R[i])>=eps) t=1;
		} 
	}while(t);
	
	// neu phuong trinh co nghiem thi thuc hien lai
    if(kt<=30){
    // sao chep lai mang ban dau
    for(i=0;i<n;i++){
	    for(j=0;j<=n;j++) A[i][j]=A2[i][j];
	    R[i]=R2[i];
	    X[i]=X2[i];
	}
	
	printf("\n");
    printf("He phuong trinh sau khi bien doi:\n");
    for(i=0;i<n;i++){
	  for(j=0;j<=n;j++)
	      printf("%8.3f",A[i][j]);
	  printf("\n");
    }
    printf("\n");
    
    printf("Qua trinh tim nghiem:\n");
    for(i=1;i<=n;i++) printf("   x%d\t",i);
    printf("\t");
    for(i=1;i<=n;i++) printf("   R%d\t",i);
    printf("\n");
	// Tim | rs | = max {|r[i]|} (i = 1->n) & tinh lai xs
	do{
		for(i=0;i<n;i++)
				printf("%8.3f",X[i]);
		printf("\t");		
		for(i=0;i<n;i++)
				printf("%8.3f",R[i]);
		printf("\n");
		t=0; // de thoat vong lap
		
		// tim max(|R[i]|)
		max=fabs(R[0]);
		k=0;
        for(i=0;i<n;i++){
        	if(max<fabs(R[i])) {
			   max=fabs(R[i]);
        	   k=i;
            }
		}	
	    X[k]+=R[k];
	    
		// tinh lai R[i] va kiem tra truong hop tiep theo
		temp=R[k];
		for(i=0;i<n;i++){
			R[i]-=A[i][k]*temp;
			if(fabs(R[i])>=eps) t=1;
		} 
	}while(t);
	//xuat nghiem
	printf("\nNghiem cua he phuong trinh la :\n");
	for(i=0;i<n;i++)
	printf("X%d=%8.3f  ",i+1,X[i]);
	printf("\n");
    }
}
//xuat ket qua ra file
void xuatfile(FILE *fileout,float X[],int n){
	fileout=fopen("output.txt","w");
	for(int i=0;i<n;i++)
	fprintf(fileout,"X%d=%8.3f  ",i+1,X[i]);
	fclose(fileout);
}
// giai bang phuong phap Gauss-Siedel
int Gauss_Siedel(float B[][Max],float X[],int n){
	printf("PHUONG PHAP GAUSS-SIEDEL:\n");
	int i,j,k,dem,lap;
	float s,A[Max][Max],A2[Max][Max],X2[Max],Y[Max],Y2[Max];
    for(i=0;i<n;i++){
	    for(j=0;j<=n;j++) A[i][j]=B[i][j];
	}
	//ham kiem tra aii co bang khong hay khong neu co thi doi dong
	aiikhackhong(A,n);
	for(i=0;i<n;i++)
	    for(j=0;j<=n;j++) A2[i][j]=A[i][j];
    printf("\nNhap vao vecto nghiem ban dau:\n  ");
    for(i=0;i<n;i++){
    	printf("x%d= ",i+1);
    	scanf("%f",&X[i]);
    	X2[i]=X[i];
	}
	// kiem tra co giai duoc khong
	dem=0;
	do{
		lap=0;dem++;
		for(i=0;i<n;i++){
			s=0;
			for(j=0;j<n;j++)
				if(i!=j) s+=A[i][j]*X[j];
			Y[i]=(A[i][n]-s)/A[i][i];
			if(fabs(Y[i]-X[i])>eps&&dem<30) lap=1;
			X[i]=Y[i];	
		}
		for (i=0; i<n; i++) X[i] = Y[i];
	}while(lap);

	if(dem==30) {
		printf("\nKhong the giai bang phuong phap nay!");
		return(0);
	}

	// neu phuong trinh co nghiem thi thuc hien lai
	else{
		// sao chep lai mang ban dau
        for(i=0;i<n;i++){
	        for(j=0;j<=n;j++) A[i][j]=A2[i][j];
	        X[i]=X2[i];
	    }

		printf("\n");
    	printf("He phuong trinh sau khi bien doi:\n");
    	for(i=0;i<n;i++){
	  		for(j=0;j<=n;j++)
	      		printf("%8.3f",A[i][j]);
	  		printf("\n");
    	}
    	printf("\n");
	    printf("Qua trinh tim nghiem:\n");
        for(i=1;i<=n;i++) printf("    x%d\t",i); 
	    printf("\n");  
	    dem=0;
	    do{
		   lap=0;dem++;
		   for (i=0;i<n;i++) printf("%8.3f", X[i]);
           printf("\n");
		   for(i=0;i<n;i++){
				s=0;
				for(j=0;j<n;j++)
					if(i!=j) s+=A[i][j]*X[j];
				Y2[i]=(A[i][n]-s)/A[i][i];
				if(fabs(Y2[i]-X[i])>eps&&dem<30) lap=1;	
			}
			for (i=0; i<n; i++) X[i] = Y2[i];
	    }while(lap);
	    //xuat nghiem
	    printf("\nNghiem cua he phuong trinh la :\n");
	    for(i=0;i<n;i++)
	        printf("X%d=%8.3f  ",i+1,X[i]);
	        printf("\n");  

	} 
}

main(){
	FILE *filein,*fileout;
	float A[Max][Max],R[Max],X[Max];
	int n,i,j,p;
	printf("***************DO AN CO SO*******************\n");
	printf("DE TAI: 906-GIAI HE PHUONG TRINH\n");
	printf("TEN THANH VIEN: NGO VAN DONG 19TCLC_DT3, TRUONG CONG VU 19TCLC_DT1\n\n");

	//mo file
	filein=fopen("input.txt","r");
	if(!filein){
		printf("Khong tim thay tep!");
		return 0 ;
	}
	
	// nhap n
	rewind(filein);
	fscanf(filein,"%d",&n);
	
	// nhap du lieu tu file va kiem tra du lieu co dung khong
	for(i=0;i<n;i++)
	    for(j=0;j<=n;j++)
		    if (!feof(filein)) 
		    	if(fscanf(filein,"%f",&A[i][j]));
				else{
				 	printf("\nSo lieu ko hop le");
	                exit(1);
				}
            else{
	            printf("\nSo lieu ko hop le");
	            exit(1);
            }
	fclose(filein);
	
	// tao menu	
	while(1){
		printf("\nMENU:\n");
	    printf("1.In ra he phuong trinh \n");
	    printf("2.Giai he phuong trinh bang phuong phap giam du\n");
		printf("3.Giai he phuong trinh bang phuong phap Gauss-Siedel\n");
	    printf("4.Xuat ket qua ra file output.txt\n");
		printf("Nhap vao so tuong ung voi chuc nang can thuc hien hoac bam so khac de thoat: ");
		scanf("%d",&p);
		switch(p){
			case 1: system("cls"); xuathe(A,n); break;
			case 2: system("cls"); giamdu(A,R,X,n); break;
            case 3: system("cls"); Gauss_Siedel(A,X,n); break;
			case 4: system("cls"); xuatfile(fileout,X,n); break;
			default: printf("\ngoodbye!"); exit(1);
		}
	}	
}
