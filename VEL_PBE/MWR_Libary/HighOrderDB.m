function HighOrderDB

%Create a DB containing the information
N=60;
HOLib_DB.GL.qx=cell(N);
HOLib_DB.GL.qw=cell(N);
HOLib_DB.GLL.qx=cell(N);
HOLib_DB.GLL.qw=cell(N);
for i=2:N
    [HOLib_DB.GL.qx{i},HOLib_DB.GL.qw{i}]=GaussLegendre(i);
    [HOLib_DB.GLL.qx{i},HOLib_DB.GLL.qw{i}]=GaussLobattoLegendre(i);
end

save('HighOrderDB','HOLib_DB')