%startar med nod l�ngst ner till v�nster, g�r motsols
n1 = [0, 0];
n2 = [0.010, 0];
n3 = [0.010, 0.010];
n4 = [0, 0.010];

node = [n1; n2; n3; n4];
%startar med elementet l�ngst ner tll h�ger

edof_nodes = [1 1 2 3; 2 3 4 1];
edof = [1 1 2 3 4 5 6; 2 5 6 7 8 1 2];

ex = [node(edof_nodes(1,2:end), 1)'; node(edof_nodes(2,2:end),1)'];
ey = [node(edof_nodes(1,2:end), 2)'; node(edof_nodes(2,2:end),2)'];

bc0 = [1 2 3 4 5 7];
bc = [bc0' bc0'*0];
dbc = -0.00004;

E = 210e9; %Är detta verkligen isotropiskt material?
t = 1;
v = 0.3;
D = hooke(2, E, v);
D = D([1 2 4], [1 2 4]);
ndof = max(max(edof));
nelm = 2;

bcMax = 0.001*3;

n_end = floor(abs(bcMax/dbc));


right_side = [3 5];
right_side_nodes = [3 5];


