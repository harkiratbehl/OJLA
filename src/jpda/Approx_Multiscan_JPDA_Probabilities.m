function Final_probabilty=Approx_Multiscan_JPDA_Probabilities(M,Obj_info,mbest)

U=size(Obj_info,1);
Final_probabilty=cell(1,U);
% Final_probabilty2=cell(1,U);
M2=sparse(M);
% h = view(biograph(M2));
% [~,C]=graphconncomp(M2);
[~,C]=graphconncomp(M2,'Directed','false');
C2=C(1:U);
NR=cell2mat(cellfun(@(x) size(x.Prob,1),Obj_info,'UniformOutput', false));

for i=unique(C2)
    ix=(C2==i);%all tracks in that component
    NR_C=NR(ix);%no of bbs(inlcuding no det) for each track
    if size(NR_C,1)==1
        TNH=NR_C;
    else
        TNH=prod(NR_C);
    end
    %         Final_probabilty(ix)=MBest_JPDA_Probabilty_Calculator(Obj_info(ix),mbest);
    %         Final_probabilty2(ix)=JPDA_Probabilty_Calculator(Obj_info(ix));
    %         erro=cell2mat(cellfun(@(x,y) any(abs(x-y)>10^-5),Final_probabilty,Final_probabilty2,'UniformOutput', false));
    %         if any(erro)
    %         error('not equal')
    %         end
    if TNH<10000||isinf(mbest)
        Final_probabilty(ix)=JPDA_Probabilty_Calculator(Obj_info(ix));
    else        
%         Final_probabilty(ix)=JPDA_Probabilty_Calculator(Obj_info(ix));
        try
            Final_probabilty(ix)=MBest_JPDA_Probabilty_Calculator(Obj_info(ix),mbest);
        catch
            Final_probabilty(ix)=JPDA_Probabilty_Calculator(Obj_info(ix));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% PROBELM COMES WHEN THIS LOOP IS EXECUTED 2 TIMES %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % which means when MBest is called
    end
    
end


