package benzenoid;

public enum CycleType {
    C5, C6, C7, UNKNOWN;

    @Override
    public String toString() {
        switch (this) {
            case C5:
                return "C5";
            case C6:
                return "C6";
            case C7:
                return "C7";
            case UNKNOWN:
            default:
                return "INCONNU";
        }
    }
}