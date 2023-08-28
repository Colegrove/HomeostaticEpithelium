package Framework.Misc;

@FunctionalInterface
public interface SweepRunFunction<T>{
    T Run(int iThread);
}
