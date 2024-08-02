def water_balance(twilt, tfield, precipitation, GAI, dates, ET0, L):
    """
    Calculate water balance for a given period.

    Parameters:
        twilt (list or array): soil wilting point values
        tfield (list or array): soil field capacity values
        precipitation (list or array): daily precipitation values
        GAI (list or array): green area index values
        dates (list of datetime.date objects): dates corresponding to each day's data
        ET0 (list or array): potential evapotranspiration values under standard conditions
        L (float): soil water content at saturation

    Returns:
        A pandas DataFrame with the calculated water balance, actual evapotranspiration,
        and date.
    """

    alpha = 0.7

    if len(precipitation) != len(GAI):
        print("Water balance function problem: GAI and precipitation have different lengths")

    length_sim = len(precipitation)

    # Initialize arrays to store results
    water = np.zeros(length_sim + 1)
    Eact = np.zeros(length_sim)
    bypass = np.zeros(length_sim)

    # Set initial water content to max field capacity
    water[0] = tfield * L

    for i in range(1, length_sim + 1):
        kc = 1.3 - 0.5 * np.exp(-0.17 * GAI[i - 1])
        ETc = ET0[i - 1] * kc
        inter = min(precipitation[i - 1], ETc, 0.2 * GAI[i - 1])  # intercepted water

        Epot = (ETc - inter)  # potential evapotranspiration
        Kr = (1 - ((0.95 * tfield[i - 1] - water[i - 1] / L) /
                   (0.95 * tfield[i - 1] - alpha * twilt[i - 1]))) ** 2

        if Kr > 1:
            Kr = 1

        Eact[i - 1] = Epot * Kr
        bypass[i - 1] = max(0, water[i - 1] - (tfield[i - 1] * L))

        water[i] = water[i - 1] + precipitation[i - 1] - Eact[i - 1] - inter - bypass[i - 1]

    # Force water content to zero if it becomes negative
    water[water < 0] = 0

    result = pd.DataFrame({'water': water[:-1], 'date': dates, 'Eact': Eact})
    result['date'] = pd.to_datetime(result['date'])
    return result