Program War_Machine;
uses Aquapoint;
Const
n= 7;
type
mas= array [1..n] of real;
mas2= array [1..n,1..n] of real;

Label Again;
var
z,pci,Tci,pcim,Tcim,pr,Tr,K,Kr,x,y,xr,yr,m,alpha,Ap,Bp,LetV,LetL,w: mas;
phi,phia,phib,phic,Bv,Av,Zv,Bl,Al,Zl,p,T,eps,sumx,sumy:real;
vA,vB,vC,vD,ppv,qqv,SSv,fiv,vA1,vB1,vC1,x1v,x2v,x3v,FFv,plbsv,y1v,y2v,y3v: real;
lA,lB,lC,lD,ppl,qql,SSl,fil,lA1,lB1,lC1,x1l,x2l,x3l,FFl,plbsl,y1l,y2l,y3l: real;
AqP:real;
Ab: mas2;
i,j:integer;
dataz,dataP,dataT,dataw:text;


function delKr(i,j: integer): integer;
begin
  result:=1;
  //if i=j then result:=0
end;

function sgn(AValue: real): integer; 
begin
  result:=0;
  if AValue>0 then result:=1;
  if AValue<0 then result:=-1;
end;

Function arch(ppvalue:real):real;
begin
  arch:= ln(ppvalue+sqrt(sqr(ppvalue)+1));
end;  

Function F(phi:real):real;
 var
 sumF:real;

 i: integer;
 begin
  sumF:=0;
   for i:= 1 to n do
    begin
     sumF:= sumF+z[i]*(K[i]-1)/(1+phi*(K[i]-1));
    End;
  F:= sumF;  
 end;
 
begin

assign (dataz, 'dataz.txt');
 reset (dataz);
for i:= 1 to n do
  read (dataz, z[i]);

assign (dataP, 'dataP.txt');  
 reset (dataP);
for i:= 1 to n do
  read (dataP, pcim[i]);

assign (dataT, 'dataT.txt');  
 reset (dataT);
for i:= 1 to n do
  read (dataT, Tcim[i]);
  
assign (dataw, 'dataw.txt');  
 reset (dataw);
for i:= 1 to n do
  read (dataw, w[i]);
 
//махинации с давл и темп  
p:=20; // атм
T:=-80;  // С

p:=p*14.6;
T:=1.8*T+491.67;  
  
//application - перевод из дурацкой метрической в нормальную,американскую.Ohh,dat feels so good
for i:= 1 to n do
begin
  pci[i]:=0.145*pcim[i];//lbs per sq_inch
  Tci[i]:=1.8*Tcim[i]+491.67;//Rankin dergees
  pr[i]:= p/pci[i];// прив давл
  Tr[i]:= T/Tci[i];// прив темп
  
 //шаг 1 - корреляция Вильсона (Вудро?Першинг?)
  K[i]:=(pci[i]/p)*exp(5.37*(1+w[i])*(1-(Tci[i]/T)));
end;  

 //еще раз перевод из метрической в нормальную,иначе паскаль сойдет с ума из-за метки
 Again:
 for i:= 1 to n do
begin
  pci[i]:=0.145*pcim[i];//lbs per sq_inch
  Tci[i]:=1.8*Tcim[i]+491.67;//Rankin dergees
  pr[i]:= p/pci[i];// прив давл
  Tr[i]:= T/Tci[i];// прив темп
end;

 //шаг 2 - Рашфор-кто-то там
  phia:= 0;
  phib:= 1;
  eps:=0.000001;
  Repeat
    phic:= (phia + phib) / 2;
    If (f(phia) * f(phic)) < 0 Then phib := phic
    Else phia := phic;
  Until abs(phib - phia) <= Eps;
   phi:= (phib+phia)/2;

 //шаг 3 - мольные доли пара и жид предварительные
 for i:= 1 to n do
 begin
  x[i]:=z[i]/(1+phi*(K[i]-1));
   y[i]:=K[i]*x[i];
 end;
 //шаг 4 - коэф.сжимаемости
  //пар
   for i:= 1 to n do
    m[i]:= 0.48+ 1.574*w[i]-0.176*w[i]*w[i];
    for i:= 1 to n do
     alpha[i]:=sqr(1+m[i]*(1-sqrt(Tr[i])));
     for i:= 1 to n do
      Ap[i]:=0.42747*alpha[i]*pr[i]/sqr(Tr[i]);
      for i:= 1 to n do
      for j:= 1 to n do
       Ab[i,j]:= sqrt(Ap[i]*Ap[j]);
       
      for i:= 1 to n do
       Bp[i]:= 0.08664*pr[i]/Tr[i];
       Av:=0;
        for i:= 1 to n do
        for j:= 1 to n do
         begin
          Av:= Av+(y[i]*y[j]*Ab[i,j]);
         end; 
       Bv:=0;
        for i:= 1 to n do
         begin
          Bv:= Bv+(y[i]*Bp[i]);
         end;
   //решение кубического уравнения .[1] неправильно считает SSv - олжен быть отр upd:все ок,правда,хз почему
    vA:= 1;
    vB:= (-1);
    vC:= Av-Bv-sqr(Bv);
    vD:=(-Av)*Bv;
     
     ppv := (3*vA*vC-vB*vB)/(3*vA*vA);
  qqv := (2*vB*vB*vB-9*vA*vB*vC+27*vA*vA*vD)/(27*vA*vA*vA);
  SSv := (4*(3*vA*vC-vB*vB)*(3*vA*vC-vB*vB)*(3*vA*vC-vB*vB)+(2*vB*vB*vB-9*vA*vB*vC+27*vA*vA*vD)*(2*vB*vB*vB-9*vA*vB*vC+27*vA*vA*vD))/(2916*vA*vA*vA*vA*vA*vA);
  if SSv<0
    then
      begin
        if qqv<0 then FFv:=Arctan(-2*Sqrt(-SSv)/qqv);
        if qqv>0 then FFv:=Arctan(-2*Sqrt(-SSv)/qqv)+Pi;
        if qqv=0 then FFv:=Pi/2;
        x1v:=2*Sqrt(-ppv/3)*Cos(FFv/3)-vB/vA/3;
        x2v:=2*Sqrt(-ppv/3)*Cos((FFv+2*Pi)/3)-vB/vA/3;
        x3v:=2*Sqrt(-ppv/3)*Cos((FFv+4*Pi)/3)-vB/vA/3;
        if qqv=0 then x3v:=-vB/vA/3;
      end;
  if SSv>0
    then
      begin
        if -qqv/2+Sqrt(SSv)>0 then y1v:=exp(ln(abs(-qqv/2+Sqrt(SSv)))/3);
        if -qqv/2+Sqrt(SSv)<0 then y1v:=-exp(ln(abs(-qqv/2+Sqrt(SSv)))/3);
        if -qqv/2+Sqrt(SSv)=0 then y1v:=0;
        if -qqv/2-Sqrt(SSV)>0 then y2v:=exp(ln(abs(-qqv/2-Sqrt(SSv)))/3);
        if -qqv/2-Sqrt(SSv)<0 then y2v:=-exp(ln(abs(-qqv/2-Sqrt(SSv)))/3);
        if -qqv/2-Sqrt(SSv)=0 then y2v:=0;
        x1v:=y1v+y2v-vB/vA/3;
        x2v:=x1v;
        x3v:= x1v;
      end;
  if SSv=0
    then
      begin
        if qqv<0 then y1v:=exp(ln(abs(-qqv/2))/3);
        if qqv>0 then y1v:=-exp(ln(abs(-qqv/2))/3);
        if qqv=0 then y1v:=0;
        x1v:=2*y1v-vB/vA/3;
        x2v:=-y1v-vB/vA/3;
        x3v:=-y1v-vB/vA/3;
      end;
       Zv:=x1v;
           if x2v>x1v then
            Zv:=x2v else
             if x3v>x1v then
              Zv:=x3v else
               Zv:=x1v;
 
  //жидкость
       Al:=0;
        for i:= 1 to n do
        for j:= 1 to n do
         begin
          Al:= Al+(x[i]*x[j]*Ab[i,j]);
         end; 
       Bl:=0;
        for i:= 1 to n do
         begin
          Bl:= Bl+(x[i]*Bp[i]);
         end;
   //решение кубического уравнения по м-ду Кардана для жидкости. [2] неправильно считает SSl - олжен быть отр upd:все ок
    lA:=1;
    lB:=-1;
    lC:=Al-Bl-sqr(Bl);
    lD:=-(Al*Bl);
     
     ppl := (3*lA*lC-lB*lB)/(3*lA*lA);
  qql := (2*lB*lB*lB-9*lA*lB*lC+27*lA*lA*lD)/(27*lA*lA*lA);
  SSl := (4*(3*lA*lC-lB*lB)*(3*lA*lC-lB*lB)*(3*lA*lC-lB*lB)+(2*lB*lB*lB-9*lA*lB*lC+27*lA*lA*lD)*(2*lB*lB*lB-9*lA*lB*lC+27*lA*lA*lD))/(2916*lA*lA*lA*lA*lA*lA);
  if SSl<0
    then
      begin
        if qql<0 then FFl:=Arctan(-2*Sqrt(-SSl)/qql);
        if qql>0 then FFl:=Arctan(-2*Sqrt(-SSl)/qql)+Pi;
        if qql=0 then FFl:=Pi/2;
        x1l:=2*Sqrt(-ppl/3)*Cos(FFl/3)-lB/lA/3;
        x2l:=2*Sqrt(-ppl/3)*Cos((FFl+2*Pi)/3)-lB/lA/3;
        x3l:=2*Sqrt(-ppl/3)*Cos((FFl+4*Pi)/3)-lB/lA/3;
        if qql=0 then x3l:=-lB/lA/3;
      end;
  if SSl>0
    then
      begin
        if -qql/2+Sqrt(SSl)>0 then y1l:=exp(ln(abs(-qql/2+Sqrt(SSl)))/3);
        if -qql/2+Sqrt(SSl)<0 then y1l:=-exp(ln(abs(-qql/2+Sqrt(SSl)))/3);
        if -qql/2+Sqrt(SSl)=0 then y1l:=0;
        if -qql/2-Sqrt(SSl)>0 then y2l:=exp(ln(abs(-qql/2-Sqrt(SSl)))/3);
        if -qql/2-Sqrt(SSl)<0 then y2l:=-exp(ln(abs(-qql/2-Sqrt(SSl)))/3);
        if -qql/2-Sqrt(SSl)=0 then y2l:=0;
        x1l:=y1l+y2l-lB/lA/3;
        x2l:=x1l;
        x3l:= x1l;
      end;
  if SSl=0
    then
      begin
        if qql<0 then y1l:=exp(ln(abs(-qql/2))/3);
        if qql>0 then y1l:=-exp(ln(abs(-qql/2))/3);
        if qql=0 then y1l:=0;
        x1l:=2*y1l-lB/lA/3;
        x2l:=-y1l-lB/lA/3;
        x3l:=-y1l-lB/lA/3;
      end;
   //Выбор коэф сжим жидкости - самый МЕНЬШИЙ из трех
   Zl:=x1l;
    if x2l<x1l then
     Zl:=x2l else
      if x3l<x1l then
       Zl:=x3l else
        Zl:=x1l;
   
 //шаг 5 - коэф летучести
  for i:= 1 to n do
   LetV[i]:= exp ((Zv-1)*Bp[i]/Bv - ln(Zv-Bv) - Av/Bv*(2*sqrt(Ap[i]/Av)-Bp[i]/Bv)*ln((Zv+Bv)/Zv)); //отриц ln(zv-bv),[1]upd:все ок, upd2:снова отрицателен
  for i:= 1 to n do
   LetL[i]:= exp ((Zl-1)*Bp[i]/Bl - ln(Zl-Bl) - Al/Bl*(2*sqrt(Ap[i]/Al)-Bp[i]/Bl)*ln((Zl+Bl)/Zl));//отриц ln(zv-bv),[2]upd:все ок,upd2:снова отрицателен
  
 //шаг 6 - новые константы равновесия и пересчет долей
  for i:= 1 to n do
   Kr[i]:=LetL[i]/LetV[i];
  for i:= 1 to n do
   xr[i]:= z[i]/(1+phi*(Kr[i]-1));
  for i:= 1 to n do   
   yr[i]:= Kr[i]*xr[i];   
 //проверка нормировка  
 for i:= 1 to n do  
   K[i]:=Kr[i];
  sumx:=0;
   for i:= 1 to n do
    sumx:= sumx+xr[i];
  sumy:=0;
   for i:= 1 to n do
    sumy:= sumy+yr[i];
   if abs(sumx-sumy)>0.001 then
    goto Again else
   writeln (sumx,' ',sumy);
   AqP:=TkCalc(z[7],P,T);
   writeln (AqP);
    for i:= 1 to n do
     begin
      writeln(i);
      writeln (Kr[i],' ',xr[i],' ',yr[i]);
     end; 

end.    