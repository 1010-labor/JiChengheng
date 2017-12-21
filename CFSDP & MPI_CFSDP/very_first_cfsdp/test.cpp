/*
#include <stdio.h>

FILE *stream;

int main(void) {

	float f;
	int c;

	FILE *fp = fopen("f:\\DATA\\data_0.dat", "r");//以只读方式打开文件。
	while ((c = fgetc(fp)) != EOF) //逐个读入字符直到文件结尾
	{
		printf("%c\n", c);
	}
	
		
		
		fseek(fp, 0L, SEEK_SET);

	
		fscanf_s(fp, "%f", &f); // C4996
	
								   // Note: fscanf is deprecated; consider using fscanf_s instead

								   // Output data read: 
	
		printf("%f\n", f);
		

		fclose(fp);
	
}*/