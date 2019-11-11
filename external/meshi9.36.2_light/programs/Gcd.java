package programs;

import java.util.Scanner;

/**
 * Prints the Greatest Common Divisor of two positive integers.
 */
public class Gcd {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int m = scanner.nextInt();
        int n = scanner.nextInt();
        int i = m;
        int j = n;

        while (i != 0) {
            int temp = j % i;
            j = i;
            i = temp;
        }

        System.out.println("Gcd = " + j);
    }
}
