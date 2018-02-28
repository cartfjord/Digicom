function [lte_gold_table] = lte_rs_gold(Ncp,Nid_cell)

% generate gold sequence

  uint32 x1=0;
  uint32 x2=0;
  
  if (Ncp==0)
    Nsymb=4;
    Ncp2=1;
  else
    Nsymb=3;
    Ncp2=0;
  end
  Nid_cellx2 = Nid_cell*2;
  for ns=0:19,
      for l=0:1,
          x2   = Ncp2 + Nid_cellx2 + (((1+Nid_cellx2)*(1+(Nsymb*l)+(7*(1+ns)))))*1024;
          fprintf('cinit (ns %d, l %d) => %d\n',ns,l,x2);
          x2_3 = bitshift(x2,-3,'uint32');
          x2_2 = bitshift(x2,-2,'uint32');
          x2_1 = bitshift(x2,-1,'uint32');
          x1   = 1+bitshift(1,31,'uint32');
          x2   = bitxor(x2,bitshift(bitxor(x2,bitxor(x2_3,bitxor(x2_2,x2_1))),31,'uint32'));
          %fprintf('n=0 : x1 %x, x2 %x\n',x1,x2);
          for n=1:49,
		  x1 = bitxor(bitshift(x1,-1,'uint32'),bitshift(x1,-4,'uint32'));
           %       fprintf('x1p : %x\n',x1);
                  x1 = bitxor(x1,bitxor(bitshift(x1,31,'uint32'),bitshift(x1,28,'uint32')));
                  x2_4 = bitshift(x2,-4,'uint32');
                  x2_3 = bitshift(x2,-3,'uint32');
                  x2_2 = bitshift(x2,-2,'uint32');
                  x2_1 = bitshift(x2,-1,'uint32');
                  x2n   = bitxor(x2_1,bitxor(x2_2,bitxor(x2_3,x2_4)));
            %      fprintf('x2p : %x\n',x2n);
                  x2 = bitxor(x2n,bitxor(bitshift(x2n,31,'uint32'),bitxor(bitshift(x2n,30,'uint32'),bitxor(bitshift(x2n,29,'uint32'),bitshift(x2n,28,'uint32'))))); 
            %      fprintf('x1 : %x, x2 : %x\n',x1,x2);
          end

          for n=0:13,
		  x1 = bitxor(bitshift(x1,-1,'uint32'),bitshift(x1,-4,'uint32'));
                  x1 = bitxor(x1,bitxor(bitshift(x1,31,'uint32'),bitshift(x1,28,'uint32')));
                  x2 = bitxor(bitshift(x2,-1,'uint32'),bitxor(bitshift(x2,-2,'uint32'),bitxor(bitshift(x2,-3,'uint32'),bitshift(x2,-4,'uint32'))));
                  x2 = bitxor(x2,bitxor(bitshift(x2,31,'uint32'),bitxor(bitshift(x2,30,'uint32'),bitxor(bitshift(x2,29,'uint32'),bitshift(x2,28,'uint32')))));
                  lte_gold_table(1+ns,1+l,1+n) = bitxor(x1,x2); 
                  fprintf('lte_gold_table[%d][%d][%d] = %x\n',ns,l,n,lte_gold_table(1+ns,1+l,1+n));
          end

      end
  end  