package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ExtendedAtomsProtein;
import meshi.parameters.ResidueType;
import meshi.util.Stat;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 14/06/12
 * Time: 14:57
 * To change this template use File | Settings | File Templates.
 */
public class Contacts {
    public static void main(String [] argv) throws IOException {
        Utils.verboseOff();
        ArrayList<Contact> contacts = getContacts();
        File dir = new File(".");
        File[] files = dir.listFiles();
        ArrayList<ContactsModel> contactsModels = new ArrayList<ContactsModel>();
        Atom atom1,atom2;
        for (File file:files) {
            boolean OK = true;
            if (file.getName().endsWith("pdb")) {
                Protein model = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(file);
                Chain chain = model.chain();
                ContactsModel contactsModel= new ContactsModel(model.name());
                for (Contact contact : contacts) {
                    Residue residue1 = chain.residueAt(contact.residueNumber1);
                    if (residue1.type() != contact.type1)
                        throw new RuntimeException("This is weird");
                    if (residue1.type() == ResidueType.GLY)
                        atom1 = residue1.ca();
                    else atom1 = residue1.cb();
                    Residue residue2 = chain.residueAt(contact.residueNumber2);
                    if (residue2.type() != contact.type2)
                        throw new RuntimeException("This is weird");
                    if (residue2.type() == ResidueType.GLY)
                        atom2 = residue2.ca();
                    else atom2 = residue2.cb();
                    double d = atom1.distanceFrom(atom2)-8;
                    if (d < 0) d = 0;
                    contactsModel.add(d);
                }
                contactsModels.add(contactsModel);
                contactsModel.setScore(Math.sqrt(contactsModel.getRms()));
                System.out.println(file.getName()+" "+contactsModel.getScore());
            }
        }
        Object[] sorted = contactsModels.toArray();
        Arrays.sort(sorted);
        System.out.println("\n          ***********");
        for (Object contactsModel: sorted) {
            System.out.println(((ContactsModel)contactsModel).name+"\t"+
                               ((ContactsModel)contactsModel).getScore());
        }
    }


    private static ArrayList<Contact> getContacts() throws IOException{
        MeshiLineReader reader = new MeshiLineReader("contacts.txt");
        String line = reader.readLine();
        String[] words = line.split(" ");
        String[] contacts = new String[(words.length-3)/2];
        for (int i = 4; i < words.length; i+=2)
             contacts[i/2-2] = words[i];
        ArrayList<Contact> out = new ArrayList<Contact>();
        for (String contactString : contacts) {
            String[] temp = contactString.split(";");
            String[] residues = temp[0].split(":");
            int l;
            int[] n = new int[2];
            ResidueType[] t = new ResidueType[2];
            for (int i = 0; i < 2; i++) {
                l = residues[i].length();
                n[i] = Integer.valueOf(residues[i].substring(0,l-1));
                t[i] = ResidueType.type(residues[i].charAt(l-1));
            }
            out.add(new Contact(n[0],t[0],n[1],t[1]));
        }
        return out;
    }
    private static class Contact {
        int residueNumber1,residueNumber2;
        ResidueType type1, type2;
        public Contact(int residueNumber1,ResidueType type1,
                      int residueNumber2,ResidueType type2){
            this.residueNumber1 = residueNumber1;
            this.residueNumber2 = residueNumber2;
            this.type1 = type1;
            this.type2 = type2;
        }
    }

    private static class ContactsModel extends Stat{
        String name;
        public ContactsModel(String name) {
            super();
            this.name = name;
        }

    }


}
