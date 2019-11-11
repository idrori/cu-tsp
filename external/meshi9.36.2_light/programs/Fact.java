package programs;

/**
 */
public class Fact {
    public static void main(String[] argv) {
        int n = 5;
        int factorial = fact(n);
        System.out.println(factorial);
        System.out.println(16 == 020);
        System.out.println(17 == 0x11);
        System.out.println(0x1A);
    }

    public static int fact(int n) {
        int ans = 1;
        while (n > 0) {
            ans = ans * n;
            n = n - 1;
        }
        return ans;
    }
}
