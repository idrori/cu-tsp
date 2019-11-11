package programs;

import java.util.Scanner;

/**


 */
public class Mod {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int m = scanner.nextInt();       // a positive integer
        int n = scanner.nextInt();       // a positive integer larger than 1
        //We will calculate n%m
        int r = m;
        while (r >= n) {
            r = r - n;
        }
        System.out.println("mode(" + m + "," + n + ") = " + r);
    }
}
