package HMMR_ATAC;

import Node.TagNode;

import java.util.*;

public class Utils {
    /**
     * Split the data by chromosome for calculation efficiency
     * @param i an ArrayList of TagNode to split
     * @return a HashMap of String and ArrayList of TagNode where the key String is the chromosome and the value ArrayList is all TagNode on that chromosome
     */
    public static HashMap<String, ArrayList<TagNode>> toMap(ArrayList<TagNode> i){
        HashMap<String,ArrayList<TagNode>> map = new HashMap<>();
        for (TagNode tagNode : i) {
            String chr = tagNode.getChrom();
            ArrayList<TagNode> temp;
            if (map.containsKey(chr)) {
                temp = map.get(chr);
            } else {
                temp = new ArrayList<>();
            }
            temp.add(tagNode);
            map.put(chr, temp);
        }

        return map;
    }

    /**
     * Split the data by chromosome for calculation efficiency, this will also clear up the original input list to save memory.
     * @param list an ArrayList of TagNode to split
     * @return a HashMap of String and ArrayList of TagNode where the key String is the chromosome and the value ArrayList is all TagNode on that chromosome
     */
    public static HashMap<String, ArrayList<TagNode>> toMapAndClear(ArrayList<TagNode> list){
        HashMap<String,ArrayList<TagNode>> map = new HashMap<>();
        Collections.reverse(list);
        ArrayList<TagNode> temp;
        for (int i = list.size()-1; i >= 0; i--){
            TagNode tagNode = list.remove(i);
            String chr = tagNode.getChrom();

            if (map.containsKey(chr)) {
                temp = map.get(chr);
            } else {
                temp = new ArrayList<>();
            }
            temp.add(tagNode);
            map.put(chr, temp);
        }
        return map;
    }
}
