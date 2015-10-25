````
2d_barrel_geometry -w 0.005 -h 0.019 -d 32 -r 0.43 -p 0.01 -n 64 -o g_test
````

````
2d_barrel_matrix -w 0.005 -h 0.019 -d 32 -r 0.43 -p 0.01 -n 64 -o m_test -e 1000000 -v --s-dl 0.04
````

````
2d_barrel_matrix -c m_test.cfg -f -o f_test m_test
````