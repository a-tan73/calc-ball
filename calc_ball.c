/* 
 * 球体の点群データから中心座標と半径を求めるプログラム
 */
#define _MAIN

// include
#include <stdio.h>
#include <math.h>
#include "mat.h"

// 定義
#define DATA_MAX 1000 // 読み込む点群データの最大数

// 構造体宣言
typedef struct
{
  float x; // x座標
  float y; // y座標
  float z; // z座標
} data_t;

int main(int argc, char *argv[])
{
  FILE *fp;

  //　変数宣言
  int i = 0;
  data_t data[DATA_MAX];
  float temp;
  double sX = 0;
  double sX2 = 0;
  double sX3 = 0;
  double sY = 0;
  double sY2 = 0;
  double sY3 = 0;
  double sZ = 0;
  double sZ2 = 0;
  double sZ3 = 0;
  double sXY = 0;
  double sYZ = 0;
  double sXZ = 0;
  double sX2Y = 0;
  double sX2Z = 0;
  double sXY2 = 0;
  double sY2Z = 0;
  double sXZ2 = 0;
  double sYZ2 = 0;
  MATRIX *Tmat, *Rmat;
  double a, b, c, d, r, z;

  // ファイル読み込み
  if ((fp = fopen(argv[1], "r")) == NULL)
  {
    printf("Don't file open.\n");
    return 0;
  }
  else
  {
    printf("open file.\n");
  }

  // データを変数に格納
  // 読み込むデータのうち、計算に必要なのは最初の３要素のみ
  while (fscanf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f",
                &data[i].x, &data[i].y, &data[i].z, &temp, &temp, &temp, &temp, &temp, &temp) == 9)
  {
    sX += data[i].x;
    sX2 += data[i].x * data[i].x;
    sX3 += data[i].x * data[i].x * data[i].x;
    sY += data[i].y;
    sY2 += data[i].y * data[i].y;
    sY3 += data[i].y * data[i].y * data[i].y;
    sZ += data[i].z;
    sZ2 += data[i].z * data[i].z;
    sZ3 += data[i].z * data[i].z * data[i].z;
    sXY += data[i].x * data[i].y;
    sXZ += data[i].x * data[i].z;
    sYZ += data[i].y * data[i].z;
    sX2Y += data[i].x * data[i].x * data[i].y;
    sX2Z += data[i].x * data[i].x * data[i].z;
    sXY2 += data[i].x * data[i].y * data[i].y;
    sY2Z += data[i].y * data[i].y * data[i].z;
    sXZ2 += data[i].x * data[i].z * data[i].z;
    sYZ2 += data[i].y * data[i].z * data[i].z;

    i++;
  }

  // 行列初期化
  double Tdat[] = {
      sX2, sXY, sXZ, sX,
      sXY, sY2, sYZ, sY,
      sXZ, sYZ, sZ2, sZ,
      sX, sY, sZ, i};
  double Rdat[] = {
      -1 * (sX3 + sXY2 + sXZ2),
      -1 * (sX2Y + sY3 + sYZ2),
      -1 * (sX2Z + sY2Z + sZ3),
      -1 * (sX2 + sY2 + sZ2)};
  Tmat = MatGenInit(4, 4, Tdat);
  Rmat = MatGenInit(1, 4, Rdat);

  // 行列計算
  Tmat = MatGenInvGJ(Tmat);
  Tmat = MatGenMul(Tmat, Rmat);

  a = Tmat->dat[0];
  b = Tmat->dat[1];
  c = Tmat->dat[2];
  d = Tmat->dat[3];

  a /= -2.0;
  b /= -2.0;
  c /= -2.0;
  r = sqrt(a * a + b * b + c * c - d);
  z = sqrt(a * a + b * b + c * c);

  // 計算結果出力
  printf("球体の中心座標(mm)　(x, y, z) = ");
  printf("(%f, %f, %f)\n", a, b, c);
  printf("半径r(mm)：%f\n", r);
  printf("レーザースキャナからの距離(mm)：%f\n", z);

  // メモリ解放
  MatFreeN(2, Tmat, Rmat);
  fclose(fp);

  return 0;
}
