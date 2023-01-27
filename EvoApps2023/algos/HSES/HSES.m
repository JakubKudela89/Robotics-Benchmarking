function res = HSES(fhd,D,Max_FEs,VarMin,VarMax)

for n=0:35
    RecordFEsFactor(n+1) = round(D^((n/10)-3.5)*Max_FEs);
end

progress = numel(RecordFEsFactor);
Run_RecordFEsFactor=RecordFEsFactor;
run_funcvals = [];
rand('state',sum(100*clock));
Xmin=VarMin;
Xmax=VarMax; bestvecs = []; bestvals = [];
for jj=1:1
    FEs=0;
    total=200;
     mu=100;
     vara=[1,-100,100,zeros(1,100),0,0];
     VRmin=repmat(Xmin,total,1);
     VRmax=repmat(Xmax,total,1);
     posinitial=VRmin+(VRmax-VRmin).*rand(total,D);
     %e=feval(fhd,posinitial',vara(:));
     [m,~] = size(posinitial);
     for i=1:m
        e(i,1)=fhd(posinitial(i,:));
     end

     FEs=FEs+total;
     weights=log(mu+1/2)-log(1:mu)';
     weights=weights/sum(weights);
     [a1,a2]=sort(e);
     bestval=a1(1); bestvals(end+1,1) = bestval;
      bestvec=posinitial(a2(1),:); bestvecs(end+1,:) = bestvec;
     pos=posinitial(a2(1:mu),:);
     meanval=mean(pos);
     stdval=std(pos);
     for k=1:total
      pos(k,:)=meanval+stdval.*randn(1,D);
     end
     for k=1:total
       for j=1:D
         if pos(k,j)>100
             pos(k,j)=meanval(j)+stdval(j).*randn;
         elseif pos(k,j)<-100
             pos(k,j)=meanval(j)+stdval(j).*randn;
         end
       end
     end
     cc1=0;
    for kk=1:100
         %e=feval(fhd,pos',vara(:));
         [m,~] = size(pos); e = zeros(m,1);
         for i=1:m
            e(i,1)=fhd(pos(i,:));
         end
         FEs=FEs+total;
             if(FEs>=Run_RecordFEsFactor(1))
                run_funcvals = [run_funcvals; bestvals(end)];
                Run_RecordFEsFactor(1)=[];
                if isempty(Run_RecordFEsFactor)
                    break;
                end
             end        
if FEs-0.01*Max_FEs<total
    processvalue(1)=a1(1);
elseif FEs-0.02*Max_FEs<total
     processvalue(2)=a1(1);
elseif FEs-0.03*Max_FEs<total
     processvalue(3)=a1(1);
     elseif FEs-0.05*Max_FEs<total
     processvalue(4)=a1(1);
     elseif FEs-0.1*Max_FEs<total
     processvalue(5)=a1(1);
     elseif FEs-0.2*Max_FEs<total
     processvalue(6)=a1(1);
     elseif FEs-0.3*Max_FEs<total
     processvalue(7)=a1(1);
     elseif FEs-0.4*Max_FEs<total
     processvalue(8)=a1(1);
elseif FEs-0.5*Max_FEs<total
     processvalue(9)=a1(1);
     elseif FEs-0.6*Max_FEs<total
     processvalue(10)=a1(1);
     elseif FEs-0.7*Max_FEs<total
     processvalue(11)=a1(1);
     elseif FEs-0.8*Max_FEs<total
     processvalue(12)=a1(1);
     elseif FEs-0.9*Max_FEs<total
     processvalue(13)=a1(1);
     elseif FEs-1*Max_FEs<total
     processvalue(14)=a1(1);
end
         
         
         [a1,a2]=sort(e);
         if a1(1)<bestval
             bestval=a1(1); bestvals(end+1,1) = bestval;
             bestvec=pos(a2(1),:); bestvecs(end+1,:) = bestvec;
         end
         newpos=pos(a2(1:mu),:);
         meanval=(newpos(:,1:D)'*weights)';
         stdval=1*std(newpos);
         FV(kk)=a1(1);
         if kk>30 
           if mod(kk,20)==0
              [aa1,aa2]=min(FV);
              if aa2<kk-20
                 cc1=1;
              end
           end
         end
         for k=1:total
             if cc1==1      %kk>300
                a=0.96*randn(1,D);
             else
                 a=randn(1,D);
             end
             pos(k,:)=meanval+stdval.*a;
         end
          for k=1:total
             for j=1:D
                if pos(k,j)>100
                   pos(k,j)=mod(pos(k,j),100);
                elseif pos(k,j)<-100
                   pos(k,j)=mod(pos(k,j),-100);
                end
             end
          end
    end
    
    previousbest=a1(1);
          if D<=30
              Times=2;
         else
              Times=1;
          end
         arfitnessbest=bestval.*ones(1,Times);
         xvalbest=repmat(bestvec',1,Times);
          N=D;
         for kkk=1:Times
         sigma=0.2;
         stopfitness=1e-8;
         if D<=30
               stopeval=Max_FEs/1;
         else
              stopeval=Max_FEs/2;
         end
         lambda=floor(3*log(N))+80;
         mu=lambda/2;
         weights=log(mu+1/2)-log(1:mu)';
         mu=floor(mu);
         weights=weights/sum(weights);
         mueff=sum(weights)^2/sum(weights.^2);
         cc=(4+mueff/N)/(N+4+2*mueff/N);
         cs=(mueff+2)/(N+mueff+5);
         c1=2/((N+1.3)^2+mueff);
         cmu=2*(mueff-2+1/mueff)/((N+2)^2+2*mueff/2);
         damps=1+2*max(0,sqrt((mueff-1)/(N+1))-1)+cs;           
          pc=zeros(N,1);
          ps=zeros(N,1);
          B=eye(N);
          DD=eye(N);
          C=B*DD*(B*DD)';
          eigenval=0;
         chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));                 %unknown
         counteval=0;
         xmean=bestvec';
         while counteval<stopeval && FEs < Max_FEs
             for k=1:lambda
                 arz(:,k)=randn(N,1);
                 arxx(:,k)=xmean+1*sigma*B*DD*arz(:,k);
              for jjj=1:N
                  if real(arxx(jjj,k))>100
                     arxx(jjj,k)=mod(real(arxx(jjj,k)),100);
                elseif real(arxx(jjj,k))<-100
                     arxx(jjj,k)=mod(real(arxx(jjj,k)),-100);
                  end
              end
                 %arfitness(k)=feval(fhd,arxx(:,k),vara(:));
                 arfitness(k)=fhd([arxx(:,k)]');
                 counteval=counteval+1;
                  FEs=FEs+1;
             end
             if(FEs>=Run_RecordFEsFactor(1))
                run_funcvals = [run_funcvals; bestvals(end)];
                Run_RecordFEsFactor(1)=[];
                if isempty(Run_RecordFEsFactor)
                    break;
                end
             end
              if FEs-0.01*Max_FEs<total
    processvalue(1)=a1(1);
elseif FEs-0.02*Max_FEs<total
     processvalue(2)=a1(1);
elseif FEs-0.03*Max_FEs<total
     processvalue(3)=a1(1);
     elseif FEs-0.05*Max_FEs<total
     processvalue(4)=a1(1);
     elseif FEs-0.1*Max_FEs<total
     processvalue(5)=a1(1);
     elseif FEs-0.2*Max_FEs<total
     processvalue(6)=a1(1);
     elseif FEs-0.3*Max_FEs<total
     processvalue(7)=a1(1);
     elseif FEs-0.4*Max_FEs<total
     processvalue(8)=a1(1);
elseif FEs-0.5*Max_FEs<total
     processvalue(9)=a1(1);
     elseif FEs-0.6*Max_FEs<total
     processvalue(10)=a1(1);
     elseif FEs-0.7*Max_FEs<total
     processvalue(11)=a1(1);
     elseif FEs-0.8*Max_FEs<total
     processvalue(12)=a1(1);
     elseif FEs-0.9*Max_FEs<total
     processvalue(13)=a1(1);
     elseif FEs-1*Max_FEs<total
     processvalue(14)=a1(1);
              end

              
                 [arfitness, arindex]=sort(arfitness);
                 xval=arxx(:,arindex(1)); 
                 
                 if abs(arfitness(1)-previousbest)<1*10^(-11)
                     break;
                 else
                     previousbest=arfitness(1);
                 end
                 
                 
                 if arfitnessbest(kkk)>arfitness(1)
                    arfitnessbest(kkk)=arfitness(1);
                    xvalbest(:,kkk)=arxx(:,arindex(1));
                 end
                 xmean=arxx(:,arindex(1:mu))*weights;
                 zmean=arz(:,arindex(1:mu))*weights;
                 ps=(1-cs)*ps+(sqrt(cs*(2-cs)*mueff))*(B*zmean);
                 hsig=norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN<1.4+2/(N+1);    
                 pc=(1-cc)*pc+hsig*sqrt(cc*(2-cc)*mueff)*(B*DD*zmean);
                 C=(1-c1-cmu)*C+c1*(pc*pc'+(1-hsig)*cc*(2-cc)*C)+cmu*(B*DD*arz(:,arindex(1:mu)))*diag(weights)*(B*DD*arz(:,arindex(1:mu)))';
                 sigma=sigma*exp((cs/damps)*(norm(ps)/chiN-1));
                 xx(counteval/lambda)=sigma;
                 if counteval-eigenval>lambda/(cmu)/N/10
                     eigenval=counteval;
                     C=triu(C)+triu(C,1)';
                     [B,DD]=eig(C);
                      DD=diag(sqrt(diag(DD)));
                 end
                 
                 if arfitness(1)==arfitness(ceil(0.7*lambda))
                    sigma=sigma*exp(0.2+cs/damps);
                    xx(counteval/lambda)=sigma;
                 end  
              %  disp([num2str(counteval) ':' num2str(arfitness(1))]);
         end
         end
         %univariate sampling 
         
         if D<=30
           total=200;
           mu=160;
         elseif D==50
           total=450;
           mu=360;
         else
              total=600;
           mu=480;
         end
        if D>=50
             if FEs<=0.3*Max_FEs
                 total=total+200;
                 mu=mu+160;
             end
         end
         weights=log(mu+1/2)-log(1:mu)';
         weights=weights/sum(weights);
         if D<=30
             ppp1=std(xvalbest');
          ppp2=sort(std(xvalbest'));
            if ppp2(1)>0.2
               dividevalue=0;
            elseif max(ppp2)<0.01
                dividevalue=1;
            else
               for dd=2:D
                   indicatorppp(dd)=(ppp2(dd)-ppp2(dd-1))/ppp2(dd-1);
               end
                   indicatorppp(1)=min(indicatorppp)-0.001;
                  [value1,value2]=sort(indicatorppp,'descend');
               for dd=1:D
                  if ppp2(value2(dd))<10
                  if ppp2(value2(dd))>0.1
                  dividevalue=ppp2(value2(dd))-0.001;
                  break;
                  end
                 elseif ppp2(value2(dd)-1)<0.01
                 dividevalue=ppp2(value2(dd))-0.001;
                break;
                  end
                  if dd==D
                       dividevalue=ppp2(value2(dd))-0.001;
                  end
                end     
            end
         else
              for kkk2=1:total/5
                  spos(kkk2,:)=xvalbest(:);
             end
         for d=1:D
            for k=1:total/5
                  spos(k,d)=xvalbest(d)-0.1*total+1*k;
            end
             [m,~] = size(spos); e = zeros(m,1);
             for i=1:m
                e(i,1)=fhd(spos(i,:));
             end
             FEs=FEs+total/5;
             bbpbbp(d)=abs(max(e)/arfitnessbest);
             for k=1:total/5
                 spos(k,d)=xvalbest(d);
             end
         end
         
          if max(bbpbbp)<3.1
             for d=1:D
                 bbpbb(d)=1;
             end
          else
              [aaa1,aaa2]=sort(bbpbbp);
              for d=1:D-1
                  diaaa1(d)=aaa1(d+1)/aaa1(d);
              end
                  [aab1,aab2]=sort(diaaa1,'descend');
              if aaa1(D/2)<=2
                for d=1:D-1
                   if aaa1(aab2(d))<1.8
                      division=aaa1(aab2(d))+0.01;
                      break;
                   end
                end
                for d=1:D
                   if bbpbbp(d)<=division
                      bbpbb(d)=1;
                   else
                       bbpbb(d)=0;
                   end
                end
              else
                  for d=1:D-1
                    if aaa1(aab2(d))<4
                      division=aaa1(aab2(d))+0.01;
                      break;
                    else division=0;
                   end
                  end
              for d=1:D
                if bbpbbp(d)<=division
                   bbpbb(d)=1;
                else
                    bbpbb(d)=0;
                end
              end
              end
          end
         end
          kk=1;
          cc2=0;
         
           VRmin=repmat(Xmin,total,1);
           VRmax=repmat(Xmax,total,1);
           pos=VRmin+(VRmax-VRmin).*rand(total,D);
           FEs=FEs+total;
           
        while FEs<Max_FEs-total
%            e1=feval(fhd,pos',vara(:));
        [m,~] = size(pos); e1 = zeros(m,1);
        for i=1:m
            e1(i)=fhd(pos(i,:));
        end
        FEs=FEs+total;
        [a1,a2]=sort(e1);
     % a1(1)
      xmin(kk)=a1(1);
        if kk==1
        [arfitnessbest,seq]=min(arfitnessbest);
        end
        if a1(1)<arfitnessbest
            xy(kk)=a1(1);
            xyvector=pos(a2(1),:);
        else
            xy(kk)=arfitnessbest;
             xyvector=xvalbest(:,seq(1));
        end
        if(FEs>=Run_RecordFEsFactor(1))
                run_funcvals = [run_funcvals; bestvals(end)];
                Run_RecordFEsFactor(1)=[];
                if isempty(Run_RecordFEsFactor)
                    break;
                end
        end 
         if FEs-0.01*Max_FEs<total
    processvalue(1)=min(xy);
elseif FEs-0.02*Max_FEs<total
     processvalue(2)=min(xy);
elseif FEs-0.03*Max_FEs<total
     processvalue(3)=a1(1);
     elseif FEs-0.05*Max_FEs<total
     processvalue(4)=min(xy);
     elseif FEs-0.1*Max_FEs<total
     processvalue(5)=min(xy);
     elseif FEs-0.2*Max_FEs<total
     processvalue(6)=min(xy);
     elseif FEs-0.3*Max_FEs<total
     processvalue(7)=min(xy);
     elseif FEs-0.4*Max_FEs<total
     processvalue(8)=min(xy);
elseif FEs-0.5*Max_FEs<total
     processvalue(9)=min(xy);
     elseif FEs-0.6*Max_FEs<total
     processvalue(10)=min(xy);
     elseif FEs-0.7*Max_FEs<total
     processvalue(11)=min(xy);
     elseif FEs-0.8*Max_FEs<total
     processvalue(12)=min(xy);
     elseif FEs-0.9*Max_FEs<total
     processvalue(13)=min(xy);
%              elseif FEs-0.99*Max_FEs<1*total
%          processvalue(14)=min(xy);
end
        
        newpos=pos(a2(1:mu),:);
        meanval=(newpos(:,1:D)'*weights)';
         stdval=1*std(newpos);
%              if max(stdval)<0.0001
%                  break;
%              end
        if kk==1
            if D>=50
               for jjj=1:D
                 if bbpbb(jjj)==0
                 stdval(jjj)=0.001;
                 meanval(jjj)=xvalbest(jjj);
                 end
               end
            else
                for jjj=1:D
                  if ppp1(jjj)<dividevalue
                      stdval(jjj)=0.001;
                 meanval(jjj)=xvalbest(jjj,seq(1));
                  end
                end
            end
        end
        kk=kk+1;
        if kk>30 
           if mod(kk,20)==0
               [aaa,bbb]=min(xmin);
           if bbb<kk-20
               cc2=1;
           else
               cc2=0;
            end
           end
        end
        for k=1:total
        if cc2==1      
               a=0.96*randn(1,D);
        else 
               a=1*randn(1,D);
        end
         pos(k,:)=meanval+stdval.*a;
        end           
         for k=1:total
             for j=1:D
                if pos(k,j)>100
                    pos(k,j)=meanval(j)+stdval(j).*randn;
               elseif pos(k,j)<-100
                    pos(k,j)=meanval(j)+stdval(j).*randn;
                end
             end
          end
%              pos=pos(1:total,:);
%              pos(1,:)=xyvector(:);
        end
     if exist('xy')
         processvalue(14)=min(xy);
         clear xy;
         clear spos;
         meanval;
     end
     finalvalue(jj)=processvalue(end);  
     process(jj,:)=processvalue;
end 
process=process';

for pr=1:14
    for pc=1:1
        if process(pr,pc)<10^(-8)
            process(pr,pc)=0;
        end
    end
end 


if finalvalue < 1e-8
    finalvalue = 0;
end
process(process< 1e-8) = 0;
res.bestval = finalvalue;
run_funcvals(36) = finalvalue;
run_funcvals(run_funcvals< 1e-8) = 0;
res.progress_values = run_funcvals;
res.progress_iters = RecordFEsFactor;
res.solution = bestvecs(end,:);
end

    

