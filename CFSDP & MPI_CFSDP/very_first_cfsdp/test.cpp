/*
#include <stdio.h>

FILE *stream;

int main(void) {

	float f;
	int c;

	FILE *fp = fopen("f:\\DATA\\data_0.dat", "r");//��ֻ����ʽ���ļ���
	while ((c = fgetc(fp)) != EOF) //��������ַ�ֱ���ļ���β
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