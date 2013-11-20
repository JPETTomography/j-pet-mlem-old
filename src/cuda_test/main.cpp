void run_kernel(char* str, int* val, int str_size, int val_size);

int main(int argc __attribute__((unused)),
         char* argv[] __attribute__((unused))) {

  char str[] = "Hello \0\0\0\0\0\0";
  int val[] = { 15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  run_kernel(str, val, sizeof(str), sizeof(val));
  return 0;
}
