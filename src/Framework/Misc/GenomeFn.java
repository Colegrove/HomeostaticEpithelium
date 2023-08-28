package Framework.Misc;

import Framework.Tools.Genome;

@FunctionalInterface
public interface GenomeFn<T extends Genome>{
    void GenomeFn(T c);
}

