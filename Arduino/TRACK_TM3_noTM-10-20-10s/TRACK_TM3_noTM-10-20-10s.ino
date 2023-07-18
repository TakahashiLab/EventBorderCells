///TM1不動　それ以外はok
#include <FlexiTimer2.h>  //Mega用タイマー割り込み
#include <avr/io.h>       //PWM出力用レジスタ設定

// Assign OUTPUT pins
#define GATE1         22
#define GATE2         23
#define GATE3         24
#define GATE4         25
#define GATE5         26
#define GATE6         27
#define GATE7         28
#define GATE8         29
#define RUN1          30  // Run treadmill 1 at speed1 FPGA 2
#define RUN2          31  // Run treadmill 2 at speed1 FPGA 1
#define DIRECTION1    32  // of Treadmill
#define DIRECTION2    33
#define FEED1         34
#define FEED2         35
#define DISCARD1      36
#define DISCARD2      37
#define FEED3           14  //changeable
#define FEED4           15  
#define DISCARD3        16
#define DISCARD4        17

// Assign INPUT pins
#define SENSOR1       38
#define SENSOR2       39
#define SENSOR3       40
#define SENSOR4       41
#define SENSOR5       42
#define SENSOR6       43
#define SENSOR7       44
#define SENSOR8       45
#define MILLSENSOR1   46                ///FPGA  青箱66
#define MILLSENSOR2   47                ///FPGA  青箱64
#define MILLSENSOR3   48                ///FPGA?  青箱62
#define MILLSENSOR4   49                ///FPGA?  青箱60
#define PORKING1        2   //changeable ///FPGA5  青箱27
#define PORKING2        3                ///FPGA6  青箱25
#define PORKING3        20               ///FPGA4  青箱19
#define PORKING4        21               ///FPGA3  青箱17
#define SWITCH        11

// state transition
#define TH  10

int prevPork=0;
bool startTask=false;
short counterTM1=0;
short counterTM2=0;
int duration=0;
short state=1;
bool ending=false;
bool endend=false;

void setup() {
  InitPinMode();
  digitalWrite(FEED1, LOW);
  digitalWrite(FEED4, LOW);
  digitalWrite(GATE4, LOW);
  digitalWrite(GATE7, HIGH);
  digitalWrite(GATE2, LOW);
  digitalWrite(GATE8, HIGH);

  pinMode(SWITCH, INPUT_PULLUP);


  pinMode(9, INPUT_PULLUP);

  //PWM波の出力用設定---------------------------------------------
  pinMode(12, OUTPUT);

  //ユーザ設定部
  unsigned int frq = 50; // 周波数←ここカメラの設定と一致させる
  float duty = 0.5; // 指定したいデューティ比
  // モード指定
  TCCR1A = 0b00000001;
  TCCR1B = 0b00010010;  //分周８
  // TOP値指定
  OCR1A = (unsigned int)(1000000 / frq);
  // Duty比指定
  OCR1B = (unsigned int)(1000000 / frq * duty);
  //-----------------------------------------------------------
  
}

void loop() {
//  if(digitalRead(9)){
//    digitalWrite(GATE1, LOW);
//  }
//  else{
//    digitalWrite(GATE1, HIGH);
//  }
//
//  
////  else if(digitalRead(8) == HIGH){
////    digitalWrite(GATE1, HIGH);
////  }
//  return;
  
  if(startTask){
    
  if (counterTM1==TH && counterTM2==TH){
    state=state+1;
    switch (state){
      case 2:
        duration=10000;
        break;
      case 3:
        duration=20000;
        break;
      case 4:
       duration=10000;
       break;
      case 5:
        ending=true;
        break;
      default:
       duration=0;
       break;
    }
    
    counterTM1=0;
    counterTM2=0;
  }
  
  if(digitalRead(PORKING3)) {
    digitalWrite(FEED3, HIGH);
    digitalWrite(FEED4, LOW);
    digitalWrite(GATE4, LOW);
    digitalWrite(GATE7, HIGH);
    digitalWrite(GATE8, HIGH);
    digitalWrite(GATE2, HIGH);
    prevPork=PORKING3;

    if(ending){
      endend=true;      
    }
  }

  
  if(digitalRead(PORKING4)) {
    digitalWrite(FEED4, HIGH);
    digitalWrite(FEED3, LOW);
    digitalWrite(GATE7, LOW);
    digitalWrite(GATE4, HIGH);
    digitalWrite(GATE2, HIGH);
    digitalWrite(GATE8, HIGH);
    prevPork=PORKING4;
     if(ending){
      endend=true;      
    }
  }

   if(digitalRead(MILLSENSOR2) && prevPork==PORKING3){
    digitalWrite(GATE2, LOW);
    digitalWrite(GATE7, LOW);
    digitalWrite(DIRECTION2, LOW);
    digitalWrite(RUN2, HIGH);
    digitalWrite(RUN1, LOW);
    FlexiTimer2::set(duration,  DelayedGateOpen2);
    FlexiTimer2::start();
    prevPork=MILLSENSOR2;
    counterTM2=counterTM2+1;
  }

   if(digitalRead(MILLSENSOR1) && prevPork==PORKING4){
    digitalWrite(GATE4, LOW);
    digitalWrite(GATE8, LOW);
    digitalWrite(DIRECTION1, LOW);
    digitalWrite(RUN1, HIGH);
    digitalWrite(RUN2, LOW);
    FlexiTimer2::set(duration,  DelayedGateOpen1);
    FlexiTimer2::start();
    prevPork=MILLSENSOR1;
    counterTM1=counterTM1+1;
  }

  if(digitalRead(SENSOR1) && prevPork==MILLSENSOR1){
    digitalWrite(GATE4, LOW);
    prevPork=MILLSENSOR1;
  }

  if(digitalRead(SENSOR5) && prevPork==MILLSENSOR2){
    digitalWrite(GATE7, LOW);
    prevPork=MILLSENSOR2;
  }
  }
  else{
    
  }
  
  //外部スイッチPWNオンオフ-------------------------------------
  if(digitalRead(9) == LOW || endend==true){
    if( TCCR1A == 0b00100001 ){
      TCCR1A = 0b00000001;///スイッチが入っているとき
    }
    else if( TCCR1A == 0b00000001 ){
      TCCR1A = 0b00100001;
    }
    delay(1000);

    if (startTask){
     digitalWrite(GATE2, LOW);
     digitalWrite(GATE4, LOW);
     ending=false;
     endend=false;
     counterTM1=0;
     counterTM2=0;
     duration=10000;
     state=1;
    }
    else{
      digitalWrite(GATE2, HIGH);
      digitalWrite(GATE4, HIGH);
    }
    startTask=!startTask;
   
  //---------------------------------------------------------
  }
  }

  void DelayedGateOpen1(){
    digitalWrite(GATE4, HIGH);
    
    digitalWrite(RUN1, LOW);
    digitalWrite(RUN2, LOW);
    FlexiTimer2::stop();
  }

  void DelayedGateOpen2(){
    digitalWrite(GATE7, HIGH);
    
    digitalWrite(RUN2, LOW);
    digitalWrite(RUN1, LOW);
    FlexiTimer2::stop();  
  }
